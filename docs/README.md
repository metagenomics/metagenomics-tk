# Meta-Omics-Toolkit

## Introduction

The Meta-Omics-Toolkit allows you to run either the full pipeline of assembly, binning and many other downstream analysis tasks or individual modules.
The toolkit can be configured by providing the module configuration via a yml file and a flag for the corresponding module or full pipeline mode.  
Options for global pipeline configuration can be viewed [here](pipeline_configuration.md).

All tools follow the same error strategy. The execution of a tool is retried three times. If the run fails the fourth time, it will be ignored.
If the execution is ignored, the toolkit will continue to run all tools that do not depend on the output of the failed tool run.
Exceptions of this handling are specified in the corresponding module section.

*Note!* Please do never place sensitive information in any of the yml configuration files since the configuration is part of the pipeline output.

## Run Full Pipeline

```BASH
nextflow run main.nf -work-dir /shared/directory/test \
	-profile PROFILE  -resume \
	-entry wPipeline -params-file example_params/fullPipeline.yml
```

where
 *  /shared/directory/test is a directory that is shared between multiple machines.
 * PROFILE can be either `standard` (local use) or `slurm` depending on which environment the pipeline should be executed.

**Note!** Metabolomics part is currently excluded from full pipeline run.


### Input

=== "Command"

    ```BASH
    -entry wPipeline -params-file example_params/fullPipeline.yml
    ```

=== "Configuration File"

    ```YAML
    ---8<--- "../example_params/fullPipeline.yml"
    ```

=== "TSV Table"

    ```TSV
    ---8<--- "../test_data/fullPipeline/reads_split.tsv"
    ```
  
    Must include the columns `SAMPLE`, `READS1` and `READS2`. `SAMPLE` must contain unique dataset identifiers
    without whitespaces or special characters. `READS1` and `READS2` are paired reads and can be HTTPS URLs, S3 links or files.

=== "Additional S3 Configuration"

    Nextflow usually stores downloaded files in the work directory. If enough scratch space is available on the worker nodes then this can be prevented by specifying
    s3 links in the input tsv file and `download` parameter in the input yaml.

    S3 TSV Links:
    ```TSV
    ---8<--- "../test_data/fullPipeline/reads_split_s3.tsv"
    ```

    Input YAML
    ```YAML
    ---8<--- "../example_params/fullPipelineQC.yml"
    ```

### Output (Overview)

In addition to the pipeline module outputs defined in the module section (Dereplication, MagAttributes, etc), the following outputs are produced. 

 * quality control (fastq) 

 * assembly (contigs)

 * binning (genomes)
 
 * read mapping (bam files)

