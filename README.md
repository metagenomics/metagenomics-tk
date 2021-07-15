
[![CI](https://github.com/pbelmann/meta-omics-toolkit/actions/workflows/workflow_modules.yml/badge.svg)](https://github.com/pbelmann/meta-omics-toolkit/actions/workflows/workflow_modules.yml)

# Meta-Omics-Toolkit

## Modules

All module configurations are the same as the full pipeline run with the sole difference that entry and param-file parameters are different.

### Run Full Pipeline

```
./nextflow run main.nf -work-dir /shared/directory/test -profile PROFILE  -resume -entry wPipeline -params-file example_params/full_pipeline.yml
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
-entry wDereplication -params-file example_params/dereplication.yml
```

#### Input

* [params-file](example_params/dereplication.yml)

* [Tsv Table](test_data/dereplication/input.tsv): Must include the columns `DATASET`, `BIN_ID`, `PATH`, `COMPLETENESS`, `CONTAMINATION`, `COVERAGE`, `N50` and `HETEROGENEITY`. 
Completeness and contamination can be used for filtering (see `params-file`). `N50`, `COVERAGE` and `HETEROGENEITY` are used for selecting the representative of every cluster.
You can set values of these columns to zero if data is not available or if you don'tt want the representative selection to be influenced by theses columns.

#### Output

The output tsv file (`final_clusters.tsv`) contains the columns `CLUSTER`, `GENOME` and `REPRESENTATIVE` where `CLUSTER` identifies a group of genomes, `GENOME` represents the path or
link of a genome and `REPRESENTATIVE` is either 0 or 1 (selected as representative). 

### MagAttributes

```
 -entry wMagAttributes -params-file example_params/magAttributes.yml 
```

#### Input

* [params-file](example_params/magAttributes.yml)

* [Tsv Table](test_data/magAttributes/input.tsv): Must include at least `DATASET` identifier and mag specific `PATH` and `BIN_ID` column.

#### Output


##### Prokka

Prokka computes `*.err`, `*.faa`, `*.ffn`, `*.fna`, `*.fsa`, `*.gbk`, `*.gff`, `*.sqn`, `*.tbl`, `*.tbl` for every bin.
Details of all files can be read on the Prokka page.
In addition it also computes a summary tsv file which adheres to the magAttributes specification.

##### GTDBTk

All GTDB files include the GTDB specific columns in addition to a `SAMPLE` column (`SAMPLE_gtdbtk.bac120.summary.tsv`, `SAMPLE_gtdbtk.ar122.summary.tsv`).
In addition, this module produces a file `SAMPLE_gtdbtk_CHUNK.tsv` that combines both files and adds a `BIN_ID` column that adheres to the magAttributes specification

##### Checkm

The Checkm output adheres to the magAttributes specification and adds a `BIN_ID` and `SAMPLE` column to the output file. 

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

## Other 

* [Developer Guidelines](docs/developer_guidelines.md)

* [Module Specification](docs/specification.md)
