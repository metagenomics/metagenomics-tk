# Global parameter settings

 * `tempdir`: Temporary directory for storing files that are used to collect intermediate files.

 * `summary`: If true a summary folder is created storing results of all samples combined.

 * `output`: Output directory for storing pipeline results. If an S3 bucket is specified with the corresponding S3 credentials (See S3 configuration section) then
   the output is written to S3.

 * `runid`: The run ID will be part of the output path and allows to distinguish between different pipeline configurations that were used for the same dataset.

 * `logDir`: A path to a directory which is used to store log files.

 * `scratch`: The scratch value can be either `false` or a path on a worker node. If a path is set, then the nextflow process in `slurm` mode is executed on the provided path.
    If the standard mode is used, then the parameter is ignored.

 * `steps`: Steps allows to specify multiple pipeline modules for running the toolkit. We distinguish between two modes. You can either run one tool of
   the pipeline or the whole pipeline with different configurations.

 * `databases`: This parameter specifies a place where files are downloaded to. If the `slurm` profile is used and databases should be downloaded, the path **should** point to a folder 
    which is not shared between the worker nodes (to reduce I/O on the shared folder resulting in a better performance).

 * `publishDirMode`: (optional) Per default results are symlinked to the chosen `output` directory. This default mode can be changed with this parameter.
    A useful mode is "copy", to copy results instead of just linking them. Other modes to choose from [here](https://www.nextflow.io/docs/latest/process.html#publishdir).  

## S3 Configuration

All modules inputs and outputs can be used in conjunction with S3.
You will have to create a configuration file that can be provided to nextflow with " -c " Parameter.
Please note that for using databases you have to provide an additional aws credentials file (see database section). 

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

## Configuration of input parameters of the full pipeline mode

### Paired End Input

The input should be a path to a tsv file containing a sample id, as well as a path to the left and right read.

Example:
```
input:
  paired:
    path: "test_data/fullPipeline/reads_split.tsv"
```

### Nanopore Input

For Nanopore data a seperate input file should be specified.

```
input:
  ont:
    path: "test_data/fullPipeline/ont.tsv"
```

### Generic SRA

The toolkit is able to fetch fastq files based on SRA run accession ids from the NCBI or from a mirror based on S3:

```
input:
  SRA:
    pattern:
      ont: ".+[^(_1|_2)].+$"
      illumina: ".+(_1|_2).+$"
    S3:
      path: test_data/SRA/samples.tsv 
      bucket: "s3://ftp.era.ebi.ac.uk" 
      prefix: "/vol1/fastq/"
      watch: false
      patternONT: ".+[^(_1|_2)].+$"
      patternIllumina: ".+(_1|_2).+$"

```

where:
  * `path` is the path to a file containing a column with `ACCESSION` as header. The `ACCESSION` column contains either SRA run or study accessions.

  * `bucket` is the S3 Bucket hosting the data.

  * `prefix` is the path to the actual SRA datasets.

  * `watch` if true, the file specified with the `path` attribute is watched and every time a new SRA run id is
     appended, the pipeline is triggered. The pipeline will never finish in this mode. Please note that watch currently only works
     if only one input type is specified (e.g "ont" or "paired" ...)

  *  `patternONT` and `patternIllumina` are patterns that are applied on the specified mirror in order to select the correct input files.

### NCBI SRA

With the following mode SRA datasets can directly be fetched from SRA.

```
input:
  SRA:
    pattern:
      ont: ".+[^(_1|_2)].+$"
      illumina: ".+(_1|_2).+$"
    NCBI:
      path: test_data/SRA/samples.tsv
```

## Database input configuration

Whenever a database field can be specified as part of the tool configuration (such as in gtdb or checkm), you are able to provide different methods to
fetch the database. In all settings, please make sure that the file has the same ending (e.g. .zip, .tar.gz) as specified in the corresponding tool section.
In addition, as database names are used to name results with which they were created, said database names should contain the respective database number or date of creation.
With this every result can be linked to one exact database version to clarify results. 
Except for the `extractedDBPath` parameter, all other input types (https, s3,...) will download the database to the folder specified in the `database` parameter.

### Extracted Database Path

If you have already downloaded and extracted the database, you can specify the path using the `extractedDBPath` parameter.
This setting is available in standard and slurm mode. In slurm mode the path can point to a db on the worker node.

Example:
```
database:
  extractedDBPath: /vol/spool/gtdb/release202
```

### HTTPS Download

The toolkit is able to download and extract the database, as long as the file ending equals the one specified in the corresponding tool section (.zip, tar.gz, tar.zst)
This setting is available in standard and slurm mode. 


Example:
```
database:
  download:
    source: 'https://openstack.cebitec.uni-bielefeld.de:8080/databases/gtdb.tar.gz'
    md5sum: 77180f6a02769e7eec6b8c22d3614d2e 
```

### Local File Path

This setting allows you to reuse an already downloaded database. 

Example:
```
database:
  download:
    source: '/vol/spool/gtdb.tar.gz'
    md5sum: 77180f6a02769e7eec6b8c22d3614d2e 
```

### S3 Download

You can specify an S3 link and configure the S3 call via the `s5cmd.params` and `s5cmd.keyfile` parameter.
The `s5cmd.params` parameter allows you to set any setting available of the [s5cmd](https://github.com/peak/s5cmd) commandline tool. 
Via the `s5cmd.kefile` parameter you can provide the path to a credentials file in case your S3 database link is not public.

In the following example the compressed file will be downloaded and extracted.

Example for publicly available compressed database:
```
database:
  download:
    source: 's3://databases/gtdb.tar.gz'
    md5sum: 77180f6a02769e7eec6b8c22d3614d2e 
    s5cmd:
      params: '--retry-count 30 --no-sign-request --no-verify-ssl --endpoint-url https://openstack.cebitec.uni-bielefeld.de:8080'
```

If your database is already extracted and available via S3, you can specify the S3 link using a wildcard as in the next example.
Example for a not publicly available uncompressed database:

```
database:
  download:
    source: 's3://databases/gtdb/*'
    md5sum: 77180f6a02769e7eec6b8c22d3614d2e 
    s5cmd:
      params: '--retry-count 30 --no-verify-ssl --endpoint-url https://openstack.cebitec.uni-bielefeld.de:8080'
      keyfile: /vol/spool/credentials
```

Your credentials file should follow the AWS credentials pattern and in case you want to use the slurm mode, be placed in a shared directory as all nodes need access:

```
[default]
aws_access_key_id=XXXXXXXX
aws_secret_access_key=XXXXXX
```

### Updating Database MD5SUMs 

The md5sum is computed over all md5sums of all files of the extracted database.
If you need to update the md5sum because you updated your database you have to download the database 
and run the following command

```
find /path/to/db -type f -exec md5sum {} \; | sort | cut -d ' ' -f 1 | md5sum | cut -d ' ' -f 1
```

### Database Download strategy

The toolkit allows to download databases on multiple nodes and tries to synchronize the download process between
multiple jobs on a node. However not all possible combinations of profiles and download types are reasonable.

| PROFILE     | Download to Shared  NFS | Download to worker scratch dir | Reuse extracted directory |
| :---------- | :--------------- | :------------------------------------ | :------------------------ |
| STANDARD    | :material-check: |  :material-close:                     | :material-check:          |
| SLURM       | :material-check-all: | :material-check:                  | :material-check:  On scratch and nfs dir |

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

The full pipeline mode is able to predict the memory consumption of some assemblers (see assembly module section).

## Fragment Recruitment for unmapped reads Configuration

Reads that could not be mapped back to a MAG can be used for fragment recruitment.
A list of genomes can be provided in the fragmentRecruitment part. 
Matched reference genomes are included in all other parts of the remaining pipeline.
Look out for their specific headers to differentiate results based on real assembled genomes and the reference genomes.
