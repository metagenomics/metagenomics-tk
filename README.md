
[![CI](https://github.com/pbelmann/meta-omics-toolkit/actions/workflows/workflow_modules.yml/badge.svg)](https://github.com/pbelmann/meta-omics-toolkit/actions/workflows/workflow_modules.yml)

# Meta-Omics-Toolkit

## Modules

All module configurations are the same as the full pipeline run with the sole difference that entry and param-file parameters are different.

### Run Full Pipeline

```
./nextflow run main.nf -work-dir /shared/directory/test -profile PROFILE  -resume -entry wPipeline -params-file example_params/full_pipeline_params.yml
```

where
 *  /shared/directory/test is a directory that is shared between multiple machines.
 * PROFILE can be either `local` or `slurm` depending on which environment the pipeline should be executed.

**Note!** Metabolomics part is currently excluded from full pipeline run.


### Run Fragment Recruitment


```
-entry wFragmentRecruitment -params-file example_params/fragmentRecruitment.yml
```

### Dereplication

```
-entry wDereplication -params-file example_params/dereplication_params.yml
```

## S3 Configuration

All modules of the pipeline can be used in conjunction with S3.
You will have to create a configuration file that can be provided to nextflow with " -c " Parameter.

```
aws {
  accessKey = xxx
  secretKey = xxx

    client {
      s_3_path_style_access = true
      connectionTimeout = 120000
      maxParallelTransfers = 28 
      maxErrorRetry = 10
      protocol = 'HTTPS'
      endpoint = 'https://openstack.cebitec.uni-bielefeld.de:8080'
      signerOverride = 'AWSS3V4SignerType'
    }
}
```

If you want to upload tool results to s3, just update the output parameter in the configuration file from `/path/to/output` to `s3://bucket/path/to/output`

# Developer Guidelines

## Testing

Tests for local use are specified in `scripts` folder. Bash scripts that start with `test_ci_` are used by github actions for continious integration tests.
Scripts for local use accept arguments for specifying local dependencies:

```
bash scripts/test_fullPipeline.sh   --steps.magAttributes.checkm.database=/vol/spool/checkm --steps.magAttributes.gtdb.database=/vol/spool/gtdb/release202
bash scripts/test_fragmentRecruitment.sh  --steps.fragmentRecruitment.frhit.genomes=test/bins/small/bin.*.fa --steps.fragmentRecruitment.frhit.samples=test/reads/small/reads.tsv 
bash scripts/test_dereplication.sh  --steps.dereplication.pasolli.input=test/bins/small/attributes.tsv
```

## Modules

Functionality is structured in modules (assembly, binning, dereplication, .etc). Each module can have multiple workflows.
Every module follows the output definition specified in the `output` section.

### Workflows

1. Worfklow names that can not be used directly and are just meant for internal use should start with an underscore.

2. At least every workflow that can be used by other external workflows should contain a short description of the functionality. 

3. Workflow names must start with `w`. 

## Process

1. Process names should start `p`

2. Processes should contain as input and output the sample id and/or the bin, contig id.


## Other

1. Magic numbers should not be used.

2. Variable, method, workflow, folder and process names should be written in camelcase.


## Output and best practice

### Motivation

* The output section is a collection of `best practices` for storing results of the `meta-omics-toolkit` output.
The definitions are motivated by the fact that the pipeline will be continuously updated and results of different pipeline
versions and modes must be differentiated.

* The idea is to run the pipeline on results of previous runs.

### Rules for dataset output

Outputs are produced by using `publish dir` directive.

```
DATASET_ID/RUN_ID/MODULE/VERSION/MODE/
```
where
   * `DATASET_ID` specifies the ID of a dataset such as the SRA dataset.
   * `RUN_ID` specifies one possible run of the full or partial pipeline. 
   * `MODULE` specifies the name of the pipeline module (e.g. binning).
   * `VERSION` specifies the module version number which follows semantic versioning (1.2.0).
   * `TOOL` specifies the name of the tool that is executed as part of the module (e.g `megahit` of the assembly module).

It is suggested that a RUN_ID output should never contain multiple versions of the same module. E.g.: 
`DATASET_ID/1/Binning/1.2.0/metabat` and `DATASET_ID/1/Binning/1.3.0/metabat`.

If a partial pipeline run (B) uses outputs of a previous run (A) (e.g. a binning tool uses the output of an assembler) and the previous run (A) alreads contains
the output of an older version of run (B), then a new RUN_ID folder must be created.

If a partial pipeline run (B) uses outputs of a previous run (A) (e.g. a binning tool uses the output of an assembler) and the previous run (A) **does not** contain
the output of an older version of run (B), then the existing RUN_ID folder must be **reused**.

### Example 1:

We assume that the following folder already exists:

```
/SRA1/1/ASSEMBLY/1.2/MEGAHIT
```


If the MODULE output does not contain a BINNING output then the existing RUN folder must be reused:

```
/SRA1/1/ASSEMBLY/1.2/MEGAHIT
/SRA1/1/BINNING/0.3/METABAT
```

### Example 2:

We assume that the following folders already exists:

```
/SRA1/1/ASSEMBLY/1.2/MEGAHIT
/SRA1/1/BINNING/0.3/METABAT
```

If the MODULE output does contain a BINNING output then a new RUN folder must be created:

```
/SRA1/1/ASSEMBLY/1.2/MEGAHIT
/SRA1/1/BINNING/0.3/METABAT
/SRA1/2/BINNING/0.4/METABAT
```

