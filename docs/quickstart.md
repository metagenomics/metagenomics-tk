# Quickstart

This quickstart allows you to run a subset of the available tools implemented by the Toolkit
on a machine to process two datasets.
This tutorial has mainly been tested on a machine with 29 GB of RAM and 14 cores running an Ubuntu operating system. 
You will need at least 250 GB of disk space. The disk were your docker images and containers are created has to be at least 17 GB.

## Requirements

1. docker: Install Docker by following the official Docker installation [instructions](https://docs.docker.com/engine/install/ubuntu/).
2. make: You can install make on Ubuntu via `sudo apt install make`.
3. java: In order to run Nextflow you need to install Java on your machine which can be achieved via `sudo apt install default-jre`.

## Preparation

Clone the metagenomics-tk repository and change the to metagenomics-tk directory.

```BASH
---8<--- "scripts/test_quickstart.sh:2:4"
```

Download the Nextflow binary that is currently supported by the Toolkit via the following command:

```BASH
---8<--- "scripts/test_quickstart.sh:6:6"
```

## Run the Toolkit

The following command will start a subset of all available modules
offered by the Toolkit. All databases will be downloaded to the database directory in your current working directory.

```BASH
---8<--- "scripts/test_quickstart.sh:8:17"
```

You can read more about the outputs, which are placed in a directory named `output`, in the corresponding [Modules](modules/introduction.md) sections.


## Further Reading

* If you want to configure these modules of want to use other datasets then paired end reads, then you can read the [configuration](configuration.md) sections.

* In case you want to directly scale out your workflow on cluster, you should continue with the [Getting Started](concept.md) section.
