This workshop demonstrates typical steps of a short-read metagenomic analysis using the Metagenomics-Toolkit and is divided into six (todo) parts.
In this part you will learn how to configure and run the Toolkit and what the output of a toolkit run looks like.

## Scope and Requirements of this Tutorial

The Metagenomics-Toolkit allows you to run either the full pipeline of assembly, binning and many other downstream analysis tasks or the individual analysis tasks.
In this tutorial you will only use the full pipeline mode. The full pipeline mode itself is structured into two parts. The first part runs the Toolkit on each
sample separately (per-sample), and the second part runs a combined downstream analysis on the output of the per-sample part, called aggregation. 
In this tutorial, you will only run the per-sample part. While there are several optimizations for running the Toolkit on a cloud-based setup, 
during this workshop you will run the Toolkit on a single machine.

### Further Requirements

* Basic Linux command line usage

* This tutorial has been tested on a machine with X cpus and X amount of RAM.  (Todo)

* Docker: Install Docker by following the official Docker installation [instructions](https://docs.docker.com/engine/install/ubuntu/).

* Java: In order to run Nextflow, you need to install Java on your machine, which can be achieved via `sudo apt install default-jre`.

* Nextflow should be installed. Please check the official Nextflow [instructions](https://www.nextflow.io/docs/latest/install.html#install-nextflow)

## Toolkit

### Execution

The Toolkit is based on Nextflow and you execute the Toolkit based on the following schema:

```BASH
NXF_VER=NEXTFLOW_VERSION nextflow run metagenomics/metagenomics-tk NEXTFLOW_OPTIONS TOOLKIT_OPTIONS 
```

* `NEXTFLOW_VERSION` is the Nextflow version that is supported by Nextflow. Every code snippet in this tutorial has a hard coded version number. 
  If you ever choose the wrong version the Toolkit will print out the versions that are supported.

* `NEXTFLOW_OPTIONS` are options that are implemented by Nextflow:

    * `-profile` determines the technology which is used to execute the Toolkit. Here we support **standard** for running the workflow on a single machine and
                 **slurm** for running the Toolkit on a cluster which uses SLURM to distribute jobs. 

	  * `-params-file` points to a configuration file that tells the Toolkit which analyses to run and which resources it should use.

	  * `-ansi-log` accepts a boolean (default: **true**) that tells Nextflow on **true** to print every update as a new line on the terminal. If **false** then Nextflow
                  prints a line for every process and updates the specific line on an update. We recommend to set **-ansi-log** to **false** because it is not possible to
                  print all possible processes on a terminal at once. 

	  * `-entry` specifies which entrypoint Nextflow should use to run the workflow. For running the full pipeline, that you will use in this workshop, you use the **wFullPipeline** entrypoint. 
               If you ever want to run seperate modules, you can check on the modules specifc page (e.g. [quality control](../../modules/qualityControl.md/#quality-control),
               [assembly](../../modules/assembly.md/#assembly))

* `TOOLKIT_OPTIONS` are options that are provided by the Toolkit. All Toolkit options are either in a configuration file which will be explained in the following or can be provided on the commandline. 

               Todo: explain how to use point notation.

### Configuration

The Toolkit uses a configuration file that specifies global parameters, the analyses that should be executed and the resources that can be used. 

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

Todo

#### Part 2: Toolkit Analyses Steps 

Most analyses, sometimes called module, that the Toolkit should execute is placed directly under the **steps** attribute.
In the example below the modules **qc** and **assembly** are directly placed under the **steps** attribute. All tools or methods
that are applied as part of the module can be seen as a property of the module. For example **MEGAHIT** is executed as part of the assembly module. 
Every level below the toolnames is for configuring the tools and methods. Every analysis is listed on the [modules page](../../modules/introduction.md) 

```YAML linenums="14" title="Example Configuration File Snippet 2"
---8<--- "default/tutorials/tutorial1/fullpipeline_assembly.yml:15:40"
```

#### Part 3: Computational Resources

The third part of a Toolkit configuration file is the **resources** attribute.
The **resources** attribute lists resource configurations where every configuration has a label and consists of a CPU count and amount of RAM.
Predefined labels are listed in the following example snippet. These labels are assigned to the processes that run the workflow specific tools. 
You can read more here about [resource parameters](../../configuration.md/#configuration-of-computational-resources-used-for-pipeline-runs).

```YAML linenums="40" title="Example Configuration File Snippet 3"
---8<--- "default/tutorials/tutorial1/fullpipeline_assembly.yml:41"
```

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

* **.command.err:** Contains standard error

* **.command.out:** Contains standard output 

* **.command.log:** Contains combined standard error and standard output

* **.command.sh:** Contains the command that was actually executed. 

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

## Workshop Data
