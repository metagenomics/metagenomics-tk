# Quickstart

This quickstart allows you to run a subset of the available tools implemented by the Toolkit
on a machine to process two datasets.
This tutorial has mainly been tested on a machine with 29 GB of RAM and 14 cores running an Ubuntu operating system. 
You will need at least 250 GB of disk space. The disk were your docker images and containers are created has to be at least 17 GB.

## Requirements

1. Docker: Install Docker by following the official Docker installation [instructions](https://docs.docker.com/engine/install/ubuntu/).
2. Java: In order to run Nextflow, you need to install Java on your machine, which can be achieved via `sudo apt install default-jre`.
3. Nextflow should be installed. Please check the official Nextflow [instructions](https://www.nextflow.io/docs/latest/install.html#install-nextflow)

## Run the Toolkit

The following command will start a subset of all available modules offered by the Toolkit. 
All databases will be downloaded to the database directory in your current working directory.

```BASH
---8<--- "scripts/test_quickstart.sh:2:14"
```

You can read more about the outputs, which are placed in a directory named `output`, in the corresponding [Modules](modules/introduction.md) sections.

## Further Reading

* If you want to configure, add or remove modules, please check the [configuration](configuration.md) section and 
check the [Getting Started](overview.md) section for an example.

* If you want to use your own datasets, then you can read the input [configuration](configuration.md#paired-end-input) sections. 
You can check the [Getting Started](overview.md) section for an example.

* In case you want to directly scale out your workflow on a cluster, you should continue with the [Getting Started](overview.md) section.
