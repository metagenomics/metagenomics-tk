# Global parameter settings

 * `tempdir`: Temporary directory for storing files that are used to collect intermediate files.

 * `summary`: If true a summary folder is created storing results of all samples combined

 * `output`: Output directory for storing pipeline results. If an S3 bucket is specified with the corresponding S3 credentials (See S3 configuration section) then
   the output is written to S3.

 * `runid`: The run ID will be part of the output path and allows to distinguish between different pipeline configurations that were used for the same dataset.

 * `logDir`: A path to a directory which is used to store log files.

 * `scratch`: The scratch value can be either `false` or a path on a worker node. If a path is set, then the nextflow process in `slurm` mode is executed on the provided path.
    If the standard mode is used, then the parameter is ignored.

 * `steps`: Steps allows to specify multiple pipeline modules for running the toolkit. We distinguish between two modes. You can either run one tool of
   the pipeline or the whole pipeline with different configurations.

## S3 Configuration

All modules of the pipeline can be used in conjunction with S3.
You will have to create a configuration file that can be provided to nextflow with " -c " Parameter.

```
aws {
  accessKey = 'xxx'
  secretKey = 'xxx'

    client {
      s_3_path_style_access = true
      connectionTimeout = 120000
      maxParallelTransfers = 28 
      maxErrorRetry = 10
      protocol = 'HTTPS'
      endpoint = 'https://openstack.cebitec.uni-bielefeld.de:8080'
      signerOverride = 'AWSS3V4SignerType'
    }
}
```

If you want to upload tool results to s3, just update the output parameter in the configuration file from `/path/to/output` to `s3://bucket/path/to/output`

If you want to use the annotation module, you also have to provide your credentials in an aws credential style.
Create a file that looks like this and fill in your credentials:

```
[default]
aws_access_key_id=ABCDEKEY
aws_secret_access_key=ABCDEKEY
```
You have to reference this file in the annotation parameter yml's s5cmd keyfiles section.  


## Configuration of input parameters of the full pipeline mode

### Generic Input

The input can be a path to a tsv file containing sample id, path to left and right read.

Example:
```
input:
  paired:
    path: "test_data/fullPipeline/reads_split.tsv"
```

### Generic SRA

The toolkit is able to fetch fastq files based on SRA run accession ids from NCBI or from a mirror based on S3:

```
input:
  SRA:
    S3:
      path: test_data/SRA/samples.tsv 
      bucket: "s3://ftp.era.ebi.ac.uk" 
      prefix: "/vol1/fastq/"
      watch: false
```

where:
  * `path` is the path to a file containing a column with `RUN_ID` as header.

  * `bucket` is the S3 Bucket hosting the data.

  * `prefix` is the path to the actual SRA datasets.

  * `watch` if true, the file specified with the `path` attribute is watched and every time a new SRA run id is
     appended, the pipeline is triggered. The pipeline will never finish in this mode.

#### NCBI SRA

With the following mode SRA datasets can directly be fetched from SRA.

```
input:
  SRA:
    NCBI:
      path: test_data/SRA/samples.tsv
```

## Optional configuration of computational resources used for pipeline runs

The toolkit uses the following machine types (flavours) for running tools. All flavours can be optionally
adjusted by modifying the cpus and memory (in GB) parameters. If for example the largest flavour is not available
in the infrastructure, `cpus` and `memory` parameters can be modified to fit the medium flavour. If larger
flavours are available, it makes especially sense to increase the `cpus` and `memory` values of the `large`
flavour to speed up for example assembly and read mapping.

Example Configuration:

```
resources:
  large:
    cpus: 28
    memory: 265
  medium:
    cpus: 14
    memory: 128
  small:
    cpus: 7
    memory: 16
  tiny:
    cpus: 1
    memory: 2
```

Additional flavours can be defined that can be used by methods that dynamically compute resources on tool error (see assembly module section).

Example:

```
resources:
  xlarge:
    cpus: 56
    memory: 512
  large:
    cpus: 28
    memory: 256
  medium:
    cpus: 14
    memory: 128
  small:
    cpus: 7
    memory: 16
  tiny:
    cpus: 1
    memory: 2
```