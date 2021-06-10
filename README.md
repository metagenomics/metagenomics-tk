
[![CI](https://github.com/pbelmann/meta-omics-toolkit/actions/workflows/workflow_modules.yml/badge.svg)](https://github.com/pbelmann/meta-omics-toolkit/actions/workflows/workflow_modules.yml)

# Meta-Omics-Toolkit

## Modules

### Run Full Pipeline

```
./nextflow run main.nf -work-dir /shared/directory/test -profile PROFILE  -resume -entry wPipeline -params-file example_params/full_pipeline_params.yml
```

where
 *  /shared/directory/test is a directory that is shared between multiple machines.
 * PROFILE can be either `local` or `slurm` depending on which environment the pipeline should be executed.

**Note!** Metabolomics part is currently excluded from full pipeline run.

### Run Fragment Recruitment

The configuration is the same as above with the sole difference that entry and param-file parameters are different:

```
-entry wFragmentRecruitment -params-file example_params/fragmentRecruitment.yml
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

## Workflows

1. Worfklow names that can not be used directly and are just meant for internal use should start with an underscore.

2. At least every workflow that can be used by other external workflows should contain a short description of the functionality. 

3. Workflow names must start with `w`. 

## Process

1. Process names should start `p`

2. Processes should contain as input and output the sample id and/or the bin, contig id.

## Other

1. Magic numbers should not be used.

2. Variable, method, workflow, folder and process names should be written in camelcase.
