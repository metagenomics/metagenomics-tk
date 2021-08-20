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
1. Create a release on Github.
2. Run `git fetch` on the master branch to get the latest tag.
2. Run `make changelog` and paste the output on the Github release section.

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
Every module follows the output definition specified in the [output specification](specification.md)  document.

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

1. Process names should start `p`

2. The input and output of processes should contain a sample and/or bin and contig id.

## Other

1. Magic numbers should not be used.

2. Variable, method, workflow, folder and process names should be written in camelcase.

