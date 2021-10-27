# Global parameter settings

 * `tempdir`: Temporary directory for storing files that are used to collect intermediate files.

 * `summary`: If true a summary folder is created storing results of all samples combined

 * `output`: Output directory for storing pipeline results. If an S3 bucket is specified with the corresponding S3 credentials (See S3 configuration section) then
   the output is written to S3.

 * `runid`: The run ID will be part of the output path and allows to distinguish between different pipeline configurations that were used for the same dataset.

 * `logDir`: A path to a directory which is used to store log files.

 * `steps`: Steps allows to specify multiple pipeline modules for running the toolkit. We distinguish between two modes. You can either run one tool of
   the pipeline or the whole pipeline with different configurations.

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
  * `path` is the path to file containing a column with `RUN_ID` as header.

  * `bucket` is the S3 Bucket hosting the data.

  * `prefix` is the path to the actual SRA datasets.

  * `watch` if true, the file specified with the `path` attribute is watched and every time a new SRA run id is
     appended, the pipeline is triggered. The pipeline will never finish in that mode.

#### NCBI SRA

With the following mode SRA datasets can directly be fetched from SRA.

```
input:
  SRA:
    NCBI:
      path: test_data/SRA/samples.tsv
```
