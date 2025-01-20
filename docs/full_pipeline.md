The full pipeline mode allows you to run the per-sample part and the aggregation part in one execution (see [schematic overview](README.md)).
In contrast to the Quickstart section, this chapter refers to the execution of the toolkit on a cluster system.
In this section, you will first run the Toolkit with a remote dataset and then learn how to replace it with your own data. 
You will then learn how to add additional analyses steps to your pipeline configuration how to add exchange tools of a module.

## Requirements

* SLURM: The Toolkit was mainly developed for cloud-based clusters using SLURM as a resource orchestrator.

* Docker: Docker must be available on all worker nodes.

* Java: In order to run Nextflow you have to install Java.

* Resources: scratch
Todo: check amount of space used???

## Part 1: Run the Toolkit with a basic config

In this section you will learn how you can run the toolkit. The input data will be downloaded automatically.
The following Toolkit analyses are performed: quality control, assembly, binning, taxonomic classification,
contamination and completeness, and annotation via Prokka.

```BASH
---8<--- "scripts/test_gettingStarted_part1.sh:2:11"
```

where

 * `NXF_HOME` points to the directory were Nextflow internal files and additional configs are stored. The default location is your home directory.
 However, it might be that your home directory is not shared among all worker nodes and is only on the master node available.  In this example
 the variable points to your current working directory (`$PWD/.nextflow`).
 * `-work-dir` points in this example to your current working directory and should point to a directory that is shared between all worker nodes.
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

Here you can see the the actual input tsv and yaml which will be provided as input.
The tsv file only describes the input data while the yaml represents the Toolkit configuration.

TODO: describe the toolkit input

=== "TSV Table"

    ```TSV
    ---8<--- "test_data/fullPipeline/reads_split.tsv"
    ```
  
    Must include the columns `SAMPLE`, `READS1` and `READS2`. `SAMPLE` must contain unique dataset identifiers
    without whitespaces or special characters. `READS1` and `READS2` are paired reads and can be HTTPS URLs, S3 links or files.

=== "Configuration File"

    ```YAML
    ---8<--- "default/fullPipeline_illumina_nanpore_getting_started_part1.yml"
    ```


### Output

The produced output can be inspected on the respective [modules page](modules/introduction.md).

## Part 2: Run the toolkit with your own data. 

Now that you have been able to run the Toolkit on your system, let's run the same configuration, but with your own data that may already be 
already on your system. We will simulate the provisioning of local files by first downloading sample paired end fastq files (size: 1,2 GB). 


```BASH
---8<--- "scripts/test_gettingStarted_part2.sh:2:4"
```

In addition, you have to tell the Toolkit the name of the sample and the path to the files that are the subject of the analysis.
For this reason, we will create a file that describes this in three columns (sample, path to left reads, path to right reads).

```BASH
---8<--- "scripts/test_gettingStarted_part2.sh:7:10"
```

Now that the files are created, you are ready to execute the Toolkit on your data but with a modified command of the first part.
The only difference is that you modify the `--input.paired.path`.

```BASH
---8<--- "scripts/test_gettingStarted_part2.sh:12:20"
```

Now you can go through the `my_data_output` folder and cheche the results.
The next section describes how to modify the analysis that was performed.

## Part 3: Exchange tools of a module 

In some cases, you may also be interested in replacing a tool from one module with another. For example, you might be interested in testing whether your selected assembler (MEGAHIT)
could work with another one like MetaSpades.

In this case you could replace the MEGAHIT part with the MetaSpades config.

```BASH
---8<--- "scripts/test_gettingStarted_part3.sh:2:17"
```

This is the MetaSpades part that is used instead of the MEGAHIT configuration:

```YAML
---8<--- "default/fullPipeline_illumina_nanpore_getting_started_part3.yml:48:51"
```

If you now compare the contigs of the two assemblers with the following command, 
you will notice that the MetaSpades assembly has a higher N50 than MEGAHIT. 

```BASH
---8<--- "scripts/test_gettingStarted_part3.sh:12:18"
```

## Part 4: Add further analysis 

In the previous sections, you have performed quality control, assembly, binning, taxonomic classification of MAGs and also basic annotation of your contigs with Prokka.
Now suppose you also want to check the presence of your genes in other databases. With the help of the Toolkit, you could create your own database with the syntax as described 
in the [wiki](database.md/#database-input-configuration). In this example, we will use the bacmet database. What you need to do here is add the following
part to the annotation section of the toolkit configuration.

The bacmet database snippet is the following

```YAML
---8<--- "default/fullPipeline_illumina_nanpore_getting_started_part4.yml:110:119"
```

By re-running the toolkit with this configuration, you will see that the previous results were cached (see next snippet) and only the annotation part is reexecuted.
Todo: mention that there is an actual default file that can be reused for all other cases.

## Further Reading

* Continue to the second part of this tutorial to learn to configure and optimize the Toolkit for your data. 

#TODO: or your infrastructure you can continue with the [configuration](configuration.md) section.

* In case you want to import the output to EMGB, please go on to the [EMGB](emgb.md) part.

* If you want to add databases or more the pre-configured ones, you can read [here](database.md) how to do this.

* You might want to adjust the [resource requirements](configuration.md/#configuration-of-computational-resources-used-for-pipeline-runs) of the Toolkit to your infrastructure.


Todo: reference other input types.
