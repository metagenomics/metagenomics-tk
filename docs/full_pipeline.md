The full pipeline mode allows you to run the per-sample part and the aggregation part in one execution (see [schematic overview](README.md)).
In contrast to the Quickstart section, this chapter refers to the execution of the Toolkit on a cluster system.
In this section, you will run the Toolkit with a dataset stored on a remote server and then learn how to replace it with your own local data. 
You will then learn how to add additional analyses steps to your pipeline configuration and how to replace the tools in a module with alternative ones.

## Requirements

1. SLURM: The Toolkit was mainly developed for cloud-based clusters using SLURM as a resource orchestrator.
2. Docker: Install Docker by following the official Docker installation [instructions](https://docs.docker.com/engine/install/ubuntu/).
3. Java: In order to run Nextflow, you need to install Java on your machine, which can be achieved via `sudo apt install default-jre`.
4. Nextflow should be installed. Please check the official Nextflow [instructions](https://www.nextflow.io/docs/latest/install.html#install-nextflow)
5. You will need at least 150 GB of scratch space on every worker node. 

## Part 1: Run the Toolkit with a basic config

In this section you will learn how to run the Toolkit. The input data will be downloaded automatically.
The following Toolkit analyses will be performed: quality control, assembly, binning, taxonomic classification,
contamination and completeness of MAGs, and gene prediction and annotation via Prokka.

!!! note "Default Configuration"
    Please note that in the following you will modify our default (best practice) configuration.  
    For most cases you don't need to modify our default configuration, you might only need to remove or add analyses.

```BASH
---8<--- "scripts/gettingStarted/test_gettingStarted_part1.sh:2:13"
```

where


 * `NXF_HOME` points to the directory where Nextflow internal files and additional configs are stored. The default location is your home directory.
 However, it might be that your home directory is not shared among all worker nodes and is only available on the master node.  In this example
 the variable points to your current working directory (`$PWD/.nextflow`).
 * `-work-dir` points in this example to your current working directory and should point to a directory that is shared between all worker nodes.
 * `-profile` defines the execution profile that should be used (local or cluster computing).
 * `-entry` is the entrypoint of the Toolkit.
 * `-params-file` sets the parameters file which defines the parameters for all tools. (see input section below)
 * `--logDir` points to a directory where your trace TSV, a timeline HTML of the executed processes and a report regarding the resource consumption of the workflow is saved.
 * `--s3SignIn` defines if any S3 login for retrieving inputs is necessary. See the [S3 configuration](configuration.md#s3-configuration) section for more information on how to configure the Toolkit for possible S3 input data.
 * `--scratch` is the directory on the worker node where all intermediate results are saved.
 * `--databases` is the directory on the worker node where all databases are saved. Already downloaded and extracted databases on a shared file system can be configured in the database setting of the corresponding [database section](database.md) in the configuration file.
 * `--output` is the output directory where all results are saved. If you want to know more about which outputs are created, then please refer to the [modules section](modules/introduction.md).
 * `--input.paired.path` is the path to a TSV file that lists the datasets that should be processed. Besides paired-end data there are also other input types. Please check the [input section](pipeline_input.md).

!!! note "Parameter override"
    Any parameters defined with a double dash are parameters that override parameters that are already specified in the YAML file.

### Input

Here you can see the actual input TSV and YAML which was used by the previous command and automatically downloaded by Nextflow. 
The TSV file only describes the input data, while the YAML file represents the Toolkit configuration.

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

You can read more about the produced output on the respective [modules page](modules/introduction.md).

## Part 2: Run the Toolkit with your own data

Now that you have been able to run the Toolkit on your system, let's use the same configuration, but with your own data that may already be 
locally available . We will simulate the provisioning of local files by first downloading sample paired-end FASTQ files (size: 1.2 GB). 
Please note that these are the same fastq files as in the previous part. The only difference to previous part is that you will download the
data beforehand. 

```BASH
---8<--- "scripts/gettingStarted/test_gettingStarted_part2.sh:2:4"
```

In addition, you have to tell the Toolkit the name of the sample and the path to the files that are the subject of the analysis.
For this reason, we will create a file that contains three columns: sample, path to left reads, and path to right reads.

```BASH
---8<--- "scripts/gettingStarted/test_gettingStarted_part2.sh:7:10"
```

Now that the files are created, you are ready to execute the Toolkit on your data but with a modified command of the first part.
The only difference is that you modify the `--input.paired.path` variable.

```BASH
---8<--- "scripts/gettingStarted/test_gettingStarted_part2.sh:12:22"
```

Now you can go through the `my_data_output` folder and check the results.
The next section describes how to modify the analysis that was performed.

## Part 3: Exchange tools of a module 

In some cases, you may also be interested in replacing a tool from one module with another. For example, you might be interested in comparing the assembler that is set as default with another one like MetaSpades.

In this case you could replace the MEGAHIT part with the MetaSpades config.

```BASH
---8<--- "scripts/gettingStarted/test_gettingStarted_part3.sh:2:13"
```

This is the MetaSpades part that was used by the previous command instead of the MEGAHIT configuration:


```YAML
---8<--- "default/fullPipeline_illumina_nanpore_getting_started_part3.yml:49:52"
```

If you now compare the contigs of the two assemblers with the following command, 
you will notice that the MetaSpades assembly has a higher N50 than the MEGAHIT one. 

```BASH
---8<--- "scripts/gettingStarted/test_gettingStarted_part3.sh:13:20"
```


## Part 4: Add further analyses 

Now, suppose you also want to check the presence of your genes in other databases. With the help of the Toolkit, you could create your own database with the syntax as described 
in the [wiki](database.md#database-input-configuration). In this example, we will use the [bacmet database](http://bacmet.biomedicine.gu.se/).
What you need to do here is to add the following part to the annotation section of the Toolkit configuration.

The bacmet database snippet is the following:

```YAML
---8<--- "default/fullPipeline_illumina_nanpore_getting_started_part4.yml:104:113"
```

By re-running the Toolkit with this configuration, you will see that the previous results were cached (see next snippet) and only the annotation part is re-executed.

```BASH
...
[a7/d8e782] Cached process > wFullPipeline:_wProcessIllumina:wShortReadBinningList:_wBinning:pMetabat (MYDATA)
[c3/537cbc] Cached process > wFullPipeline:_wProcessIllumina:wShortReadBinningList:_wBinning:pCovermContigsCoverage (Sample: MYDATA)
[bb/fcd1aa] Cached process > wFullPipeline:_wProcessIllumina:wShortReadBinningList:_wBinning:pCovermGenomeCoverage (Sample: MYDATA)
[8c/08e157] Cached process > wFullPipeline:wMagAttributesList:_wMagAttributes:pGtdbtk (Sample: MYDATA)
[4b/286c49] Cached process > wFullPipeline:wMagAttributesList:_wMagAttributes:pCheckM2 (Sample: MYDATA)
...
```

## Further Reading

* Continue to the aggregation part of the Getting Started [tutorial](aggregation.md) to learn how to aggregate data of multiple samples.

* You can check in our [configuration](configuration.md) section for how to further adapt the Toolkit to your infrastructure.

* In case you want to import the output to EMGB, please visit the [EMGB](modules/export.md) section.

* If you want to add your own sequence databases, you can read [here](database.md) how to do this.

* You might want to adjust the [resource requirements](configuration.md#configuration-of-computational-resources-used-for-pipeline-runs) of the Toolkit to your infrastructure.
