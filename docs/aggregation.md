There are two ways to execute the Toolkit. You can either run all steps in one execution or you run first the per-sample analysis
(e.g. assembly, binning, annotation, etc.) and afterwards combine the results (e.g. via dereplication and co-occurrence) in a second run.
The second option allows you to process multiple samples via independent Toolkit executions on different infrastructures and combine all
results afterwards.

## Requirements

1. SLURM: The Toolkit was mainly developed for cloud-based clusters using SLURM as a resource orchestrator.
2. Docker: Install Docker by following the official Docker installation [instructions](https://docs.docker.com/engine/install/ubuntu/).
3. Java: In order to run Nextflow you need to install Java on your machine which can be achieved via `sudo apt install default-jre`.
4. Nextflow should be installed. Please check the official Nextflow [instructions](https://www.nextflow.io/docs/latest/install.html#install-nextflow)
5. This tutorial assumes that you have already executed the Toolkit as described in the [full pipeline section](full_pipeline.md).

## Run the Toolkit

```BASH
---8<--- "scripts/gettingStarted/test_gettingStarted_aggregation.sh:1:13"
```

where

 * `-work-dir` points to a directory that is shared between multiple machines.
 * `-profile` defines the execution profile that should be used.
 * `-entry` is the entry point of the aggregation workflow.
 * `-params-file` sets the parameters file which defines the parameters for all tools. (see input section below)
 * `--logDir` points to a directory where your trace tsv, a timeline html of the executed processes and a report regarding the resource consumption of the workflow is saved.
 * `--s3SignIn` defines if any S3 login for retrieving inputs is necessary. See the [S3 configuration](configuration.md/#s3-configuration) section for more information on how to configure the Toolkit for possible S3 input data.
 * `--scratch` is the directory on the worker node where all intermediate results are saved.
 * `--databases` is the directory on the worker node where all databases are saved. Already downloaded databases on a shared file system can be configured in the database setting of the corresponding [database section](database.md) in the configuration file.
 * `--output` is the output directory where all results are saved. If you want to know more about which outputs are created, then please refer to the [modules section](modules/introduction.md).
 * `--input` points to the output directory of the per-sample workflow.

!!! note "Parameter override"
    Any parameters defined with a double dash are parameters that override parameters that are already specified in the yaml file.

### Input

=== "Configuration File"

    ```YAML
    ---8<--- "default/fullPipelineAggregate.yml"
    ```

### Output

The meaning of the produced output can be inspected on the respective [module page](modules/introduction.md).
You can check for the results in the `AGGREGATED` folder:

For example you could check for the number of species clusters created through dereplication:

```BASH
cat  my_data_spades_output/AGGREGATED/1/dereplication/*/bottomUpClustering/clusters/clusters.tsv
```

!!! note "Parameter override"
    Please note that the dereplication method produces more meaningful results when more then one sample is
    provided as input.

## Further Reading

* Pipeline Configuration: If you want to configure and optimize the Toolkit for your data or your infrastructure then you can continue with the [configuration](configuration.md) section.

* In case you want to import the output to EMGB, please go on to the [EMGB](modules/export.md) configuration section. Please keep in mind that for EMGB only the per-sample part is necessary.

* You might want to adjust the [resource requirements](configuration.md/#configuration-of-computational-resources-used-for-pipeline-runs) of the Toolkit to your infrastructure.
