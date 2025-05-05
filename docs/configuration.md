## Global parameter settings

 * `tempdir`: Temporary directory for storing files that are used to collect intermediate files.

 * `output`: Output directory for storing pipeline results. If an S3 bucket is specified with the corresponding S3 credentials (See S3 configuration section) then
   the output is written to S3.

 * `runid`: The run ID will be part of the output path and allows to distinguish between different pipeline configurations that were used for the same dataset.

 * `logDir`: A path to a directory which is used to store log files.

 * `scratch`: The scratch value can be either `false` or a path on a worker node. If a path is set, then the nextflow process in `slurm` mode is executed on the provided path.
    If the standard mode is used, then the parameter is ignored.

 * `steps`: Steps allows to specify multiple pipeline modules for running the toolkit. We distinguish between two modes. You can either run one tool of
   the pipeline or the whole pipeline with different configurations.

 * `databases`: This parameter specifies a place where files are downloaded to. If the `slurm` profile is used and databases should be downloaded, the path **should** point to a folder 
    which is not shared between the worker nodes (to reduce I/O on the shared folder resulting in a better performance). If this parameter is provided, the toolkit will create the specified
    directory. If all your databases have already been extracted beforehand, you can simply omit this parameter.

 * `publishDirMode`: (optional) Per default results are symlinked to the chosen `output` directory. This default mode can be changed with this parameter.
    A useful mode is "copy", to copy results instead of just linking them. Other modes to choose from [here](https://www.nextflow.io/docs/latest/process.html#publishdir).  

 * `skipVersionCheck`: The toolkit is regurarly tested against a set of Nextflow versions. Setting the `--skipVersionCheck` allows you to use the toolkit with Nextflow versions
   that were not tested.

 * `s3SignIn`: If your input data (not the databases) is not publicly accessible via S3, then you will have to set the `s3SignIn` parameter to `true`.

## S3 Configuration

All module inputs and outputs can be used in conjunction with S3.
If you want to set a custom S3 configuration setting (i.e. custom S3 endpoint), you will have to modify the aws client parameters 
with " -c ".

Example:
```
---8<--- "test_data/assets/aws.config"

```

In addition you will have to set a Nextflow Secret with the following keys:

```
nextflow secrets set S3_ACCESS xxxxxxxxx
nextflow secrets set S3_SECRET xxxxxxxxx
```

`S3_ACCESS` corresponds to the aws S3 access key id and `S3_SECRET` is the aws S3 secret key.
If your input data (not the databases) is publicly available then you have to set `s3SignIn:` to `false` in your config file.
Please note that for using databases you have to set additional credentials ([see database section](database.md/#s3-download)). 

## Configuration of Computational Resources used for Pipeline Runs

The toolkit uses the following machine types (flavors) for running tools. All flavors can be optionally
adjusted by modifying the cpus and memory (in GB) parameters. If for example the largest flavor is not available
in the infrastructure, `cpus` and `memory` parameters can be modified to fit the highmemMedium flavor. If larger
flavors are available, it makes especially sense to increase the `cpus` and `memory` values of the `large`
flavor to speed up for example assembly and read mapping.

Example Configuration:

```
resources:
  highmemLarge:
    cpus: 28
    memory: 230
  highmemMedium:
    cpus: 14
    memory: 113
  large:
    cpus: 28
    memory: 58
  medium:
    cpus: 14
    memory: 29
  small:
    cpus: 7
    memory: 14
  tiny:
    cpus: 1
    memory: 1
```

Additional flavors can be defined that can be used by methods that dynamically compute resources on tool error (for example the [assembly module](modules/assembly.md)).

Example:

```
resources:
  xlarge:
    cpus: 56
    memory: 512
  highmemLarge:
    cpus: 28
    memory: 230
  highmemMedium:
    cpus: 14
    memory: 113
  large:
    cpus: 28
    memory: 58
  medium:
    cpus: 14
    memory: 29
  small:
    cpus: 7
    memory: 14
  tiny:
    cpus: 1
    memory: 1
```

The full pipeline mode is able to predict the memory consumption of some assemblers (see [assembly module](modules/assembly.md)). The prediction
will also consider additional flavors which have been added to the resources section in the configuration file.

## Apptainer (Experimental)

Apptainer allows containers to be executed without relying on a daemon process, making it a suitable software solution for HPC setups.
Currently, the Toolkit only supports a limited set of processes that can be run using Apptainer, and even those are only tested using the local execution on one machine (`standard` profile).

Apptainer was installed for testing on an Ubuntu machine by following the instructions on the official [website](https://apptainer.org/docs/admin/latest/installation.html#install-ubuntu-packages).
Here, the `apptainer-suid` package was installed.

The following command executes the Toolkit on a local machine using Apptainer. The configuration corresponds to the one on the Quickstart page,
except that Smetana, which is part of the metabolomics module, is disabled.

```BASH linenums="1" title="Quickstart command for using Apptainer"
    ---8<--- "scripts/test_quickstart_apptainer.sh:3:13"

