This workshop demonstrates typical steps of a short-read metagenomic analysis using the Metagenomics-Toolkit and is divided into eight parts.
In this part you will learn how to configure and run the Toolkit and what the output of a Toolkit run looks like.

## Tutorial Scope and Requirements

The Metagenomics-Toolkit allows you to run either the full pipeline of assembly, binning and many other downstream analysis tasks or the individual analyses seperately. 
In this tutorial you will only use the *full pipeline* mode. The *full pipeline* mode itself is structured into two parts. The first part runs the Toolkit on each
sample separately (*per-sample*), and the second part runs a combined downstream analysis on the output of the *per-sample* part, called *aggregation*. 
In this tutorial, you will only run the *per-sample* part. While there are several optimizations for running the Toolkit on a cloud-based setup, 
during this workshop you will run the Toolkit on a single machine.

### Requirements

!!! warning "Course Participants"

    In case you are a course participant, it is very likely that a machine has been prepared for you and you can ignore this section.
    If in doubt, ask the course organizers. 

* Basic Linux command line usage

* This tutorial has been tested on a machine with 28 cpus and 64 GBs of RAM with Ubuntu installed on it.

* Docker: Install Docker by following the official Docker installation [instructions](https://docs.docker.com/engine/install/ubuntu/).

* Java: In order to run Nextflow, you need to install Java on your machine, which can be achieved via `sudo apt install default-jre`.

* Nextflow should be installed. Please check the official Nextflow [instructions](https://www.nextflow.io/docs/latest/install.html#install-nextflow).

### Preparations

For this tutorial, we need some preparations on our machine.
**Please note that the following instructions can be customized for your specific machine.** 

#### Link to your volume

We are using the location `~/mgcourse/` as directory, where all the analysis is stored. We assume that there is a larger storage volume mounted at `/vol/mgcourse`. 
If your preferred volume is at another location, change the `/vol/mgcourse` to a different path accordingly. If you don't have a another volume with more storage available
you can also simply leave out the following command and create the directory with `mkdir ~/mgcourse`.

!!! warning "Course Participants"

    In case you are a course participant, a volume has been prepared for you, so please run the following command. 


!!! question "Task 1"
    We will link the volume to `~/mgcourse`, so we can use that link for the rest of this tutorial:
    ```BASH
    sudo ln -s /vol/mgcourse/ ~/mgcourse
    ```

#### Database directory link

In addition, we need a symlink for the directory where the databases should be stored.

!!! question "Task 2"
    Create a `databases` and `/vol/scratch` directory and create a link to the database directory:
    ```BASH
    sudo mkdir ~/mgcourse/databases
    sudo mkdir -p /vol/scratch
    sudo ln -s ~/mgcourse/databases/ /vol/scratch/databases
    ```
 

## Metagenomics-Toolkit Introduction

### Execution

The Toolkit is based on Nextflow and you can execute the Toolkit based on the following commandline schema:

```BASH
NXF_VER=NEXTFLOW_VERSION nextflow run metagenomics/metagenomics-tk NEXTFLOW_OPTIONS TOOLKIT_OPTIONS 
```

* `NEXTFLOW_VERSION` is the Nextflow version that is supported by Nextflow. Every code snippet in this tutorial has a hard coded version number. 
  If you ever choose the wrong version, then the Toolkit will print out the versions that are supported.

* `NEXTFLOW_OPTIONS` are options that are implemented by Nextflow:

    * `-profile` determines the technology which is used to execute the Toolkit. Here we support **standard** for running the workflow on a single machine and
                 **slurm** for running the Toolkit on a cluster which uses SLURM to distribute jobs. 

	  * `-params-file` points to a configuration file that tells the Toolkit which analyses to run and which resources it should use. An example configuration file will be explained in the next section.

	  * `-ansi-log` accepts a boolean (default: **true**) that tells Nextflow on **true** to print every update as a new line on the terminal. If **false** then Nextflow
                  prints a line for every process and updates the specific line on an update. We recommend to set **-ansi-log** to **false** because it is not possible to
                  print all possible processes on a terminal at once when running the Toolkit.

	  * `-entry` specifies which entrypoint Nextflow should use to run the workflow. For running the *full pipeline*, that you will use in this workshop, you use the **wFullPipeline** entrypoint. 
               If you ever want to run seperate modules, you can check on the modules specifc page (e.g. [assembly](../../modules/assembly.md/#assembly)).

* `TOOLKIT_OPTIONS` are options that are provided by the Toolkit. All Toolkit options are either in a configuration file or can be provided on the commandline which will be explained in the following section. 


!!! question "Task 3"

    Open the Metagenomics-Toolkit wiki on a second browser tab by a click on this
    [link](https://metagenomics.github.io/metagenomics-tk/latest/){:target="_blank"}.
    Imagine you need only to run the quality-control part seperately. Can you tell the name of the **entrypoint**? 
    Use the wiki page you have opened on another tab to answer the question.

    ??? Solution
        If you go to the [quality control](../../modules/qualityControl.md/#quality-control) part, then you will find
        the **wShortReadQualityControl** entrypoint for short reads and the **wOntQualityControl** entrypoint for long
        reads.

### Configuration

The Toolkit uses a YAML configuration file that specifies global parameters, the analyses that will be executed and the computational resources that can be used. 

The configuration file is divided into three parts:

#### Part 1: Global Workflow Parameters

The following snippet shows parameters that effect the whole execution of the 
workflow. All parameters are explained in a dedicated Toolkit wiki [section](../../configuration.md). 

```YAML linenums="1" title="Example Configuration File Snippet 1"
---8<--- "default/tutorials/tutorial1/fullpipeline_assembly.yml:1:14"
```

!!! info "Computational Resources"
    Please note that computational resources are also global parameters and will be handled in the third part of this configuration section. 

##### Input Field

The input field (line 3, snippet 1) specifies the type of input data to process (Nanopore, Illumina, data hosted on SRA or a mirror)
and you can find a dedicated wiki section [here](../../pipeline_input.md). Regardless of which input type
is used, the user must provide a file containing a list of datasets to be processed.
The list can be a list of remote or local files and in the case of SRA, a list of SRA run IDs.

Since you will work with short read data in this tutorial, your input file looks like this:

```BASH linenums="1" title="Input File"
---8<--- "test_data/tutorials/tutorial1/reads.tsv"
```

The first columns (SAMPLE) specifies the name of the dataset. The second (READS1) and third (READS2) column specify the files
containing the forward and reverse reads.

#### Part 2: Toolkit Analyses Steps 

Analyses, or sometimes called modules, that the Toolkit executes are placed directly under the **steps** attribute in the configuration file.
In the example below, the modules **qc** and **assembly** are placed directly under the **steps** attribute. Any tools or methods
that are used as part of the module can be considered as a property of the module. For example **MEGAHIT** is executed as part of the assembly module.
The level below the tool names is for configuring the tools and methods. Each analysis is listed on the [modules page](../../modules/introduction.md). 

```YAML linenums="14" title="Example Configuration File Snippet 2"
---8<--- "default/tutorials/tutorial1/fullpipeline_assembly.yml:15:40"
```

#### Part 3: Computational Resources

The third part of a Toolkit configuration file is the **resources** attribute.
The **resources** attribute lists computational resource configurations, where each configuration has a label and consists of the number of CPUs and amount of RAM assigned to it.
Predefined labels are listed in the following example snippet. These labels are assigned to the processes that run the workflow specific tools.
You can read more  about resource parameters [here](../../configuration.md/#configuration-of-computational-resources-used-for-pipeline-runs).

```YAML linenums="40" title="Example Configuration File Snippet 3"
---8<--- "default/tutorials/tutorial1/fullpipeline_assembly.yml:41"
```

!!! question "Task 4"

    One of the first checks before running the Toolkit is to adjust the resource labels to the resources of your machine.
    You can run `nproc` to get the number of CPUs and `free -h --mega` to get the amount of RAM (row name: Mem, column name: total) available on your machine. 
    Is there enough RAM on your machine to run the Toolkit? 

    ??? Solution
        If you are using a machine as described in the Requirements section, then yes, there is enough RAM available for the workflow. 


#### Configuration File vs. Commandline Parameters

All parameters defined in the YAML configuration file can also be supplied as command-line arguments. To do this, prefix each parameter with a double dash (--). 
If a parameter is nested within the hierarchy of the YAML file, represent it as a command-line argument by connecting each level of the hierarchy using a dot (.).

For example, consider the CPU count of the *highmemLarge* resource label in the previous snippet.
The corresponding command-line argument would be `--resources.highmemLarge.cpus`.

!!! question "Task 5"

    Let`s say you want to specify a path to a different input tsv file (see **Example Configuration File Snippet 1**) that contains a different set of input datasets. 
    How would you specify the parameter on the commandline?

    ??? Solution
        `--input.paired.path`


### Output

The Toolkit output fulfills the following schema:

```BASH
SAMPLE_NAME/RUN_ID/MODULE/MODULE_VERSION/TOOL
```

* **RUN_ID:** The run ID will be part of the output path and allows to distinguish between different pipeline configurations that were used for the same dataset.

* **MODULE** is the analysis that is executed by the Toolkit (e.g. qc, assembly, etc.). 

* **MODULE_VERSION** is the version number of the module.

* **TOOL** is the tool or method that is executed by the Toolkit.

Below you can see an example output structure.
Every output folder includes four log files:

* **.command.err:** Contains the standard error.

* **.command.out:** Contains the standard output. 

* **.command.log:** Contains the combined standard error and standard output.

* **.command.sh:** Contains the command that was executed. 

```BASH linenums="1" title="Example Output Directory"
output/
└── sample
    └── 1
        └── qc
            └── 0.3.0
                ├── fastp
                │   ├── .command.err
                │   ├── .command.log
                │   ├── .command.out
                │   ├── .command.sh
                │   ├── sample_fastp.json
                │   ├── sample_fastp_summary_after.tsv
                │   ├── sample_fastp_summary_before.tsv
                │   ├── sample_interleaved.qc.fq.gz
                │   ├── sample_report.html
                │   ├── sample_unpaired.qc.fq.gz
                │   └── sample_unpaired_summary.tsv
                ├── filterHuman
                │   ├── .command.err
                │   ├── .command.log
                │   ├── .command.out
                │   ├── .command.sh
                │   ├── sample_interleaved.filtered.fq.gz
                │   ├── sample_interleaved.removed.fq.gz
                │   ├── sample_interleaved_summary_after.tsv
                │   ├── sample_interleaved_summary_before.tsv
                │   ├── sample_unpaired.filtered.fq.gz
                │   ├── sample_unpaired_summary_after.tsv
                │   └── sample_unpaired_summary_before.tsv
                ├── kmc
                │   ├── .command.err
                │   ├── .command.log
                │   ├── .command.out
                │   ├── .command.sh
                │   ├── sample.13.histo.tsv
                │   ├── sample.13.kmc.json
                │   ├── sample.21.histo.tsv
                │   ├── sample.21.kmc.json
                │   ├── sample.71.histo.tsv
                │   └── sample.71.kmc.json
                └── nonpareil
                    ├── .command.err
                    ├── .command.log
                    ├── .command.out
                    ├── .command.sh
                    ├── sample.npa
                    ├── sample.npc
                    ├── sample.npl
                    ├── sample_nonpareil_curves.pdf
                    └── sample_nonpareil_index.tsv

```

  


