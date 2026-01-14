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
      skipDB: false
```

where:
  * `path` is the path to a file containing a column with `ACCESSION` as header. The `ACCESSION` column contains either SRA run or study accessions.

  * `bucket` is the S3 Bucket hosting the data.

  * `prefix` is the path to the actual SRA datasets.

  * `watch` if true, the file specified with the `path` attribute is watched and every time a new SRA run id is
     appended, the pipeline is triggered. The pipeline will never finish in this mode. Please note that watch currently only works
     if only one input type is specified (e.g "ont" or "paired" ...)

  *  `patternONT` and `patternIllumina` are patterns that are applied on the specified mirror in order to select the correct input files.

  *  `skipDB` is a flag for skipping SRA run IDs from the local SRA database. This flag is only used for debugging.

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

### Configurtion of input parameters of the aggregation mode

```
input:
  perSampleOutput: "output"
  selectedSamples: "test_data/fullPipeline/filter.tsv" 
```

where:
  * `perSampleOutput` is the output folder of the per sample run

  * `selectedSamples` is an optional parameter that allows you to select specific samples of interest.
  The output of these samples is located in the `perSampleOutput` directory. This option is useful when not all the samples in an output directory are to be used as input for modules such as Read Mapping or Cooccurrence.
