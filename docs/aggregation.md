There are two ways to execute the toolkit. You can either run all steps in one execution or you run first the per sample analysis
(e.g. assembly, binning, annotation, etc.) and afterwards you combine the results (e.g. dereplication, co-occurrence) in a second run.
The second option allows you to process multiple samples via independent toolkit executions on different infrastructures and combine all
results afterwards.

You could run first the wFullPipeline mode as described in the [full pipeline section](full_pipeline.md) but without dereplication,
read mapping and co-occurrence modules and afterwards run the the aggregation as described below:

## Requirements

* SLURM: The Toolkit was mainly developed for cloud-based clusters using SLURM as a resource orchestrator.

* Docker: Docker must be available on all worker nodes.

* Java: In order to run Nextflow you have install Java.

* Nextflow version: Please check in the nextflow.config file the supported Nextflow versions.

* Resources: scratch

## Preparation

First checkout the Github repository in a directory which is shared by all worker nodes:

## Run the Toolkit

```BASH
---8<--- "scripts/test_gettingStarted_aggregation.sh:1:11"
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
 * `--input` points to the output directory of the per-sample workflow.

!!! note ""
    The output directory is the same directory that is used as input.

!!! note "Parameter override"
    Any parameters defined with a double dash are parameters that override parameters that are already specified in the yaml file.


### Input

=== "Configuration File"

    ```YAML
    ---8<--- "default/fullPipelineAggregate.yml"
    ```

### Output

The produced output can be inspected on on the [modules page](modules/introduction.md).

## Further Reading


* Pipeline Configuration: If you want to configure and optimize the Toolkit for your data or your infrastructure you can continue with the [configuration](configuration.md) section.

* In case you want to import the output to EMGB, please go on to the [EMGB](emgb.md) part. Please keep in mind that for EMGB only the per-sample part is necessary.

* If you want to add databases or more the pre-configured ones, you can read [here](database.md) how to do this.

* You might want to adjust the [resource requirements](configuration.md/#configuration-of-computational-resources-used-for-pipeline-runs) of the Toolkit to your infrastructure.
