The full pipeline mode allows you to run the per-sample part and the aggregation part in one execution (see [schematic overview](README.md)).
In contrast to the quickstart section, this chapter refers to the execution of the toolkit on a cluster system.

## Requirements

* SLURM: The Toolkit was mainly developed for cloud-based clusters using SLURM as a resource orchestrator.

* Docker: Docker must be available on all worker nodes.

* Java: In order to run Nextflow you have install Java.

* Nextflow version: Please check in the nextflow.config file the supported Nextflow versions.

* Resources: scratch

## Preparation

First checkout the Github repository in a directory which is shared by all worker nodes:

```BASH
git clone git@github.com:metagenomics/metagenomics-tk.git
cd metagenomics-tk
```

## Run the Toolkit

In the following case the Nextflow binary is placed on a share working directory:

```BASH
./nextflow run main.nf -work-dir /vol/spool/work \
    -profile slurm \
    -entry wFullPipeline \
    -params-file default/fullPipeline_illumina_nanpore_getting_started.yml  \
    --s3SignIn false \
    --scratch /vol/scratch \
    --databases /vol/scratch/databases \
    --output output \
    --input.paired.path test_data/fullPipeline/reads_split.tsv
```

where

 * `-work-dir` points to a directory that is shared between multiple machines.
 * `-profile` defines the execution profile that should be used.
 * `-entry` is the entrypoint of the Toolkit.
 * `-params-file` sets the parameters file which defines the parameters for all tools. (see input section below)
 * `--s3SignIn` defines if any S3 login for retrieving inputs is necessary. See the [S3 configuration](configuration.md/#s3-configuration) section for more information on how to configure the Toolkit for possible S3 input data.
 * `--scratch` is the directory on the worker node where all intermediate results are saved.
 * `--databases` is the directory on the worker node where all databases are saved. Already downloaded databases on a shared file system can be configured in the database setting of the corresponding [database section](database.md) in the configuration file.
 * `--output` is the output directory where all results are saved. If you want to know more about which outputs are created, then please refer to the [modules section](modules/introduction.md).
 * `--input.paired.path` is the path to a tsv file that lists the datasets that should be processed. If you want to provide other parameters, please check the [input section](pipeline_input.md).

!!! note "Parameter override"
    Any parameters defined with a double dash are parameters that override parameters that are already specified in the yaml file.


### Input

=== "TSV Table"

    ```TSV
    ---8<--- "test_data/fullPipeline/reads_split.tsv"
    ```
  
    Must include the columns `SAMPLE`, `READS1` and `READS2`. `SAMPLE` must contain unique dataset identifiers
    without whitespaces or special characters. `READS1` and `READS2` are paired reads and can be HTTPS URLs, S3 links or files.

=== "Configuration File"

    ```YAML
    ---8<--- "default/fullPipeline_illumina_nanpore.yml"
    ```

=== "Additional S3 Configuration"

    Nextflow usually stores downloaded files in the work directory. If enough scratch disc space is available on the worker nodes then this can be prevented by specifying
    s3 links in the input tsv file and `download` parameter in the qc part of the input yaml.

    !!! note ""
     
        **Currently the S3 links must point to files that are publicly available.**

    S3 TSV Links:
    ```TSV
    ---8<--- "test_data/fullPipeline/reads_split_s3.tsv"
    ```

    Example Input YAML with the `download` parameter in the qc part.
    ```YAML
    ---8<--- "example_params/fullPipelineQC.yml"
    ```


### Output

The produced output can be inspected on on the [modules page](modules/introduction.md).

## Further Reading

* Pipeline Configuration: If you want to configure and optimize the Toolkit for your data or your infrastructure you can continue with the [configuration](configuration.md) section.

* In case you want to import the output to EMGB, please go on to the [EMGB](emgb.md) part.

* If you want to add databases or more the pre-configured ones, you can read [here](database.md) how to do this.

* You might want to adjust the [resource requirements](configuration.md/#configuration-of-computational-resources-used-for-pipeline-runs) of the Toolkit to your infrastructure.
