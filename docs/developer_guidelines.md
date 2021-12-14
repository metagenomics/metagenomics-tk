# Developer Guidelines

## Commit and Release Guidelines

We are using [git-chglog](https://github.com/git-chglog/git-chglog) to automatically generate a changelog for the latest released based on our commit messages.
Commit messages should follow the following format:

```
feat(scope): feature added in the scope
```

Example:

```
feat(assembly): megahit added
```

`feat` can be replaced by one of the formats specified in the options sections of the [config](.chglog/config.yml) file.
Scope can for example represent a module, a configuration or a specific document.

A new release should be made the following way: 
1. Update pipeline version in the nextflow manifest `nextflow.config`.
2. Create a release on Github.
3. Run `git fetch` on the master branch to get the latest tag.
4. Run `make changelog` and paste the output on the Github release section.

## Testing

Tests for local use are specified in the `scripts` folder. Bash scripts that start with `test_ci_` are used by github actions for continious integration tests.
Scripts for local use accept arguments for specifying local dependencies:

Examples:
```
bash scripts/test_fullPipeline.sh   --steps.magAttributes.checkm.database=/vol/spool/checkm --steps.magAttributes.gtdb.database=/vol/spool/gtdb/release202
bash scripts/test_fragmentRecruitment.sh  --steps.fragmentRecruitment.frhit.genomes=test/bins/small/bin.*.fa --steps.fragmentRecruitment.frhit.samples=test/reads/small/reads.tsv 
bash scripts/test_dereplication.sh  --steps.dereplication.pasolli.input=test/bins/small/attributes.tsv
bash scripts/test_magAttributes.sh  --steps.magAttributes.input=test/bins/small/attributes.tsv
```

## Modules

Functionality is structured in modules (assembly, binning, dereplication, .etc). Each module can have multiple workflows.
Every module follows the output definition specified in the [output specification](specification.md)  document. The name and the version of the
module is specified in the `modules` section of the `nextflow.config` file.

### Workflows

1. Worfklow names that can not be used directly and are just meant for internal use should start with an underscore.

2. At least every workflow that can be used by other external workflows should contain a short description of the functionality. 

3. Workflow names must start with `w`. 

## Versioning

All modules are versioned according to [semantic versioning](https://semver.org/). The version number is incorporated in the output directory (see [output specification](specification.md)) 
for easier parsing of the output directory. In the following we give examples when to increment which part of of the version identifier:

Given a version number MAJOR.MINOR.PATCH, increment the:

  * MAJOR version when you make incompatible changes, as for example modifying the output structure. A script that was build to parse the output structure must be adapted then.
  * MINOR version when you add functionality in a backwards compatible manner. One example is adding an additional tool to the module. 
  * PATCH version when you make backwards compatible bug fixes. This is necessary when you for example increment the docker container version number that fixes a bug or increases the
    speed of the tool.

## Process

1. Process names should start `p`.

2. The input and output of processes should contain a sample and/or bin and contig id.

3. Pocesses should publish `.command.sh`, `.command.out`, `.command.log` and `.command.err` files but never `.command.run`.
In cases where processes process different data but publish it to the same folder these files would be overwritten on every run.
For example when Prokka publishes log files of every genome to the same sample directory.
For that reason these files need to be renamed, so that their names include a unique id (e.g. bin id). 
Please output those files to channel with the following entries and connect this channel to the pDumpLogs process that you can import
from the utils module:

```
include { pDumpLogs } from '../utils/processes'

...

tuple env(FILE_ID), val("${output}"), val(params.LOG_LEVELS.INFO), file(".command.sh"), \
        file(".command.out"), file(".command.err"), file(".command.log"), emit: logs
```

Examples can be viewed in the Checkm and Prokka process.

## Logs

Log files should be stored in the user provided `logDir` directory.

### Log Level

Every configuration file must have a `logLevel` attribute that can have the following values:

```
ALL = 0  All logs are published
INFO = 1 Just necessary logs are published
```

These values can be used in the publish dir directive to enable or disable the output of logs.


```
   publishDir params.output, saveAs: { filename -> getOutput(params.runid, "pasolli/mash/sketch", filename) }, \
        pattern: "{**.out,**.err, **.sh, **.log}", enabled: params.logLevel <= params.LOG_LEVELS.ALL
```

Furthermore the `params.LOG_LEVELS.*` parameters can be used inside of a process to enable or disable intermediate results for debugging purposes.
In case the log is the send to the pDumpLogs process (see Process section), you can specify the log level as part of the tuple:

```
tuple env(FILE_ID), val("${output}"), val(params.LOG_LEVELS.INFO), file(".command.sh"), \
        file(".command.out"), file(".command.err"), file(".command.log"), emit: logs
```

4. Custom error strategies that do not follow the strategy defined in nextflow.config, should be documented (see Megahit example).

## Databases

<<<<<<< HEAD
If the same databases is downloaded during runtime by multiple processes, it takes up an unnecessary ammount of disc space.
One idea is to always use the same place to store these databases. This place should be described in `params.databases`.
If other processes try to use this databses they can look at `params.databases` on this current machine. 
If it is present it can be used, if not it should be downloaded. Through this procedure only one copy of each databases 
is used, which is space-saving.   
=======
If the same database is downloaded during runtime by multiple processes, it takes up an unnecessary ammount of disc space.
One idea is too always use the same place to store these databases. This place should be described in `params.databases`.
If other processes try to use this databases they can look at `params.databases` on the current machine. 
If it is present it can be used, if not it should be downloaded. Through this procedure only one copy of each databases is used,
which is space-saving.   
>>>>>>> master

## Configuration

Every process should be configurable by providing a parameters string to the tool in the process.
Every module should use the following specification in the configuration file:

```
steps:
  moduleName:
    parameter: 42
    processName:
      additionalParams: " --super-flag "
```

Additional params can have a string value (like the example above) that is provided to the tool:

```
pProcess {

   ...

  shell:
  """
  supertool !{params.steps.moduleName.processName.parameter}  !{params.steps.moduleName.processName.additionalParams}
  """
}

```

The value of the `additionalParams` key can also be a map if multiple tools are used in the same process:

``
steps:
  moduleName:
    parameter: 42
    processName:
      additionalParams:
         toolNameA: " -c 84  "
         toolNameB: " --super-flag "
```

`parameter` fields can hold hardcoded parameters that hold a defined value like a number that should not be a string.
One use case of those parameters is that they can be reused for multiple tools.

Example:

```
pProcess {

   ...

  shell:
  """
  toolNameA --super-specific-number-flag !{params.steps.moduleName.parameter}
  toolNameB --similar-flag-to-toolA !{params.steps.moduleName.parameter} 
  """
}

```

### Internal Configuration

The `_wConfigurePipeline` workflow in the main.nf file should be used for setting 
pipeline parameters that are need for fullfilling the user provided configuration.

Example:
Lets assume the user enables the plasmid module. In that case it is mandatory that 
the assembler produces a fastg file independend of the user provided settings of the assembler.
The `_wConfigurePipeline` method sets in that case the fastg parameter of any assembler to `true`.

## Tools

### concurrentDownload.sh

This script allows to synchronize the download of a database between multiple jobs and should be executed the following way.
Before a database is downloaded, the script checks either the MD5SUM or the database version against a user specified parameter.

```
flock LOCK_FILE concurrentDownload.sh --output=DATABASE \
           --httpsCommand=COMMAND \
           --localCommand=COMMAND \
           --s3Command=COMMAND \
           --link=LINK \
           --expectedMD5SUM=USER_VERIFIED_DATABASE_MD5SUM
```

where
  * `LOCK_FILE` is a file that is used for locking. Processes will check if the file is currently locked before trying to download anything.
    This file should ideally placed in the `params.database` directory of the specific tool (e.g. !{params.databases}/rgi).
  
  * `DATABASE` is the directory that is used for placing the specific database.

  * `COMMAND` is the command used to download and extract the database. (e.g. "wget -O data $DOWNLOAD_LINK && tar -xvf data ./card.json && rm data" for the `--httpsCommand` flag)

  * `USER_VERIFIED_DATABASE_MD5SUM` is the MD5SUM of the *extracted* database that the user should test manually before executing the pipeline.

  * `LINK` is the link that will be used to test if is S3, HTTPS or a local path.

## Other

1. Magic numbers should not be used.

2. Variable, method, workflow, folder and process names should be written in camelcase.
