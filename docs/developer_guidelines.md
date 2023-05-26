# Guidelines

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

`feat` can be replaced by one of the formats specified in the options sections of the config file (see example below).
Scope can for example represent a module, a configuration or a specific document.

A new release should be made the following way: 

1. Update pipeline version in the nextflow manifest `nextflow.config`.

2. Create a release on Github.

3. Run `git fetch` on the master branch to get the latest tag.

4. Run `make changelog` and paste the output on the Github release section.

```YAML
---8<--- "../.chglog/config.yml"
```

### Versioning

Following [semantic versioning](https://semver.org/), we define the configuration input file and the output folder structure as our public `API`.
Changes to the version numbers reflect updates to the config or the output folders and files.
The toolkit consists of many modules that can be used in different combinations and
because of this flexibility, we had to come up with a detailed versioning system. We version each module separately, as well as the pipeline itself.

Module version numbers are updated when a module-specific input parameter is updated or the output folder or file structure is changed.
All module version numbers can be retrieved by running the toolkit with the `wGetModuleVersion` entrypoint and should be reported on the release page. 

The module version number is incorporated in the output directory (see [output specification](pipeline_specification.md)) 
for easier parsing of the output directory. In the following we give examples when to increment which part of of the version identifier:

Given a version number MAJOR.MINOR.PATCH, increment the:

  * MAJOR version when you make incompatible changes, as for example modifying the output structure. A script that was build to parse the output structure must be adapted then.
  * MINOR version when you add functionality in a backwards compatible manner. One example is adding an additional tool to the module. 
  * PATCH version when you make backwards compatible bug fixes. This is necessary when you for example increment the docker container version number that fixes a bug or increases the
    speed of the tool.

The pipeline specific version number defined in the mainifest part of the nextflow.config should be changed
if either any module specific version number is incremented or any module-independent parameter (e.g. `tempdir`) or output structure is changed. 

## Testing

Tests for local use are specified in the `scripts` folder. These scripts are also used as part of the continuous integration tests.
If you want to run these scripts locally, you will have to override the paths to the databases you have downloaded:

Examples:
```BASH
bash scripts/test_fullPipeline.sh  " --steps.magAttributes.checkm.database=/vol/spool/checkm --steps.magAttributes.gtdb.database=/vol/spool/gtdb/release202 "
bash scripts/test_fragmentRecruitment.sh  " --steps.fragmentRecruitment.frhit.genomes=test/bins/small/bin.*.fa --steps.fragmentRecruitment.frhit.samples=test/reads/small/reads.tsv "
bash scripts/test_dereplication.sh "  --steps.dereplication.pasolli.input=test/bins/small/attributes.tsv "
bash scripts/test_magAttributes.sh "  --steps.magAttributes.input=test/bins/small/attributes.tsv "
```

### Nextflow Versions

The toolkit is tested against the lowest and highest Nextflow version number specified in VERSIONS.txt.

## Modules

Functionality is structured in modules (assembly, binning, dereplication, .etc). Each module can have multiple workflows.
Every module follows the output definition specified in the [output specification](pipeline_specification.md)  document. The name and the version of the
module is specified in the `modules` section of the `nextflow.config` file.

### Workflows

1. Worfklow names that can not be used directly and are just meant for internal use should start with an underscore.

2. At least every workflow that can be used by other external workflows should contain a short description of the functionality. 

3. Workflow names must start with `w`. 

## Process

Process names should start `p`. The in- and output of processes should contain a sample and/or a bin and contig id.
Custom error strategies that do not follow the strategy defined in nextflow.config, should be documented (see Megahit example).

### Processes should publish process specific files

Processes should publish `.command.sh`, `.command.out`, `.command.log` and `.command.err` files but never `.command.run`.
In cases where processes process different data but publish it to the same folder these files would be overwritten on every run.
For example when Prokka publishes log files of every genome to the same sample directory.
For that reason these files need to be renamed, so that their names include a unique id (e.g. bin id). 
Please output those files to channel with the following entries and connect this channel to the pDumpLogs process that you can import
from the utils module:

```JAVA
include { pDumpLogs } from '../utils/processes'

...

tuple env(FILE_ID), val("${output}"), val(params.LOG_LEVELS.INFO), file(".command.sh"), \
        file(".command.out"), file(".command.err"), file(".command.log"), emit: logs
```

Examples can be viewed in the Checkm and Prokka process.

#### Logs

Log files should be stored in the user provided `logDir` directory.

##### Log Level

Every configuration file must have a `logLevel` attribute that can have the following values:

```
ALL = 0  All logs are published
INFO = 1 Just necessary logs are published
```

These values can be used in the publish dir directive to enable or disable the output of logs.


```JAVA
   publishDir params.output, mode: "${params.publishDirMode}", saveAs: { filename -> getOutput(params.runid, "pasolli/mash/sketch", filename) }, \
        pattern: "{**.out,**.err, **.sh, **.log}", enabled: params.logLevel <= params.LOG_LEVELS.ALL
```

Furthermore the `params.LOG_LEVELS.*` parameters can be used inside of a process to enable or disable intermediate results for debugging purposes.
In cases where the log is send to the pDumpLogs process (see Process section), you can specify the log level as part of the tuple:

```
tuple env(FILE_ID), val("${output}"), val(params.LOG_LEVELS.INFO), file(".command.sh"), \
        file(".command.out"), file(".command.err"), file(".command.log"), emit: logs
```

### Time Limit

Every process must define a time limit which will never be reached on "normal" execution. This limit is only useful for errors in the execution environment  
which could lead to an endless execution of the process.

You can use the setTimeLimit helper method to add a user configurable time limit.

Example:

```
time Utils.setTimeLimit(params.steps.qc.fastp, params.modules.qc.process.fastp.defaults, params.mediumDefault)
```


## Databases

If the same database is downloaded during runtime by multiple processes, it takes up an unnecessary ammount of disc space.
One idea is too always use the same place to store these databases. This place should be described in `params.databases`.
If other processes try to use this databases they can look at `params.databases` on the current machine. 
If it is present it can be used, if not it should be downloaded. Through this procedure only one copy of each databases is used,
which is space-saving. Links to the actual database should contain the database version number or the date of download.

## Configuration

Every process should be configurable by providing a parameters string to the tool in the process.
Every module should use the following specification in the configuration file:

```YAML
steps:
  moduleName:
    parameter: 42
    processName:
      additionalParams: " --super-flag "
      timeLimit: "AUTO"
```

Please check the process chapter regarding possible values for the time limit attribute.
Additional params can have a string value (like the example above) that is provided to the tool:

```JAVA
pProcess {

   ...

  shell:
  """
  supertool !{params.steps.moduleName.processName.parameter}  !{params.steps.moduleName.processName.additionalParams}
  """
}

```

The value of the `additionalParams` key can also be a map if multiple tools are used in the same process:

```YAML
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
In that case the fastg parameter of any assembler will be set to `true` by the `_wConfigurePipeline` method.

## Toolkit Docker Images

Dockerfiles of Docker images that are build by toolkit developers can be found in the `docker` directory.
The name of the directory (i.e.: `toolkit-python-env` in `docker/toolkit-python-env`) is used for the docker image name. All images belong to the
[metagenomics](https://quay.io/organization/metagenomics) quay.io organisation which is owned by the Computational Metagenomics group in Bielefeld.
A docker repository in the `metagenomics` orginsation must be created by the organisation owner, before the actual image can be build.
The version of the image specified in the `VERSION` file (i.e. `docker/toolkit-python-env/VERSION`) is used for the image tag (`metagenomics/toolkit-python-env:VERSION`).
An image build is only triggered if the version in the VERSION file is updated on the dev or master branch.

## Wiki

For building the documentation we are using [mkdocs](https://www.mkdocs.org/) in combination with [mkdocs-material](https://squidfunk.github.io/mkdocs-material/)
and a [plugin](https://timvink.github.io/mkdocs-print-site-plugin/index.html) for building static single page html files.
The wiki HTML files are uploaded to S3 storage on pull request merge events in the master and dev branch (see Makefile commands using `make help`).

You can work on these html files locally by running `make dev_wiki`. But please note that by build the static html file for upload, the navigation might change.
You can view the final html file by building the html file (see Makefile `make help`). 

## Utils

We do not want to duplicate code and thats why we should store methods in the lib/Utils.groovy file. The Utils class can be used in any module. 

## Database Download

This section explains how a developer is able to implement the database download strategy as explained in the [user documentation](pipeline_configuration.md#database-input-configuration). 
Example implementations can be found in the gtdb, checkm or rgi scripts.

The first step is to check if the user provides an already extracted database: 

```BASH
DB_PATH=""
if [ -z "!{EXTRACTED_DB}" ]
then
   # Fetch user parameter for not extracted db path and run flock (see next section)
   DB_PATH="not extracted"
else
  # Set variable to extracted db path
fi
```

Since the download is not directly handled by nextflow and paths to the files need to be downloaded, any file or directory must be
mounted first to the container. For this reason you have to add the `setDockerMount` function with the database config as input to 
the `containerOptions` parameter:

```BASH
containerOptions " other container options " + setDockerMount(params.steps?.magAttributes?.checkm?.database)
```

### Filesystem locks

Multiple jobs of the same process (e.g. GTDB) are able to synchronize the download of a database by using filesystem locks.
The download is handled by the `concurrentDownload.sh` script and should be executed the following way:

```BASH
flock LOCK_FILE concurrentDownload.sh --output=DATABASE \
           --httpsCommand=COMMAND \
           --localCommand=COMMAND \
           --s3FileCommand=COMMAND \
           --s3DirectoryCommand=COMMAND \
           --s5cmdAdditionalParams=S5CMD_PARAMS \
           --link=LINK \
           --expectedMD5SUM=USER_VERIFIED_DATABASE_MD5SUM
```

where
  * `LOCK_FILE` is a file that is used for locking. Processes will check if the file is currently locked before trying to download anything.
    This file should ideally placed in the `params.database` directory of the specific tool (e.g. !{params.databases}/rgi).
  
  * `DATABASE` is the directory that is used for placing the specific database.

  * `COMMAND` is the command used to download and extract the database and to remove it afterwards. 
    (e.g. "wget -O data.tar.gz $DOWNLOAD_LINK && tar -xvf data.tar.gz ./card.json && rm data.tar.gz" for the `--httpsCommand` flag)

  * `USER_VERIFIED_DATABASE_MD5SUM` is the MD5SUM of the *extracted* database that the user should test manually before executing the pipeline.

  * `S5CMD_PARAMS` allows you to set s5cmd specific parameters. For more information check the s5cmd documentation. 

  * `LINK` is the link that will be used to test if the file is accessible by S3, HTTPS or is available via a local path.

  * `USER_VERIFIED_DATABASE_MD5SUM` Before a database is downloaded, the script checks the MD5SUM of an already downloaded database against a user specified one.
     If it does not equal, the script will download the database again.

### Tests

You can test your tool against different database inputs by using the `make runDatabaseTest` command. You will have to specify multiple databases 
that are accessible via https, S3, local path etc. Please check github actions file for how to run these tests.

## Polished Variables

Sometimes user input variables must be polished before they can used in our code.
Thats why the nextflow config adds a namespace to the params namespace called `polished`.
For example the params.databases variable must end with a slash in order to be used as part of a docker mount.
Thats why there is a variable `params.polished.databases` that should be used instead.  

## Other

1. Magic numbers should not be used.

2. Variable, method, workflow, folder and process names should be written in camelcase.

