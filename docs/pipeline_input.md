## Configuration of input parameters of the full pipeline mode

### Paired End

#### Sample Sheet

The input can be a path to a tsv file containing a sample id, as well as a path to the left and right read.

Example:
```YAML
input:
  paired:
    sheet: "test_data/fullPipeline/reads_split.tsv"
```

The sample sheet must have the columns `SAMPLE`, `READS1` and `READS2`.
The `SAMPLE` column must be unique and the `READS1` and `READS2` columns point to the repective left and right read files.
`READS1` and `READS2` can be local paths, URLs or S3 paths. In case S3 is used, additional [configuration](pipeline_input.md#s3-Configuration) is necessary.

Example:
```
---8<--- "test_data/fullPipeline/reads_split.tsv"

```

#### Command line

For a small number of samples it is sometimes easier to provide them directly on the command line instead of creating a sample sheet. 
In that case you can provide the input as follows:

```BASH
	  --input.paired.r1 read1.fq.gz \
	  --input.paired.r2 read2.fq.gz \
	  --input.paired.names test1
```

where

  * The `--input.paired.r1` and `--input.paired.r2` parameters point to the left and right reads, respectively.
    The left and right reads of the respective sample must be provided in the same order and the number of files must be the same.

  * `--input.paired.names` are the sample names. If multiple samples are provided they must be enclosed in double ticks (e.g. "test1 test2").
    The number of sample names must match the number of files provided for `--input.paired.r1` and `--input.paired.r2`.

The `--input.paired.r1` and `--input.paired.r2` can point to the same type of resources (URL, S3, etc.) as the `READS1` and `READS2` columns (see "Sample Sheet" section).
See [Quickstart](quickstart.md#run-the-toolkit) section for a working example.

### Nanopore

#### Sample Sheet

For Nanopore data a seperate sample sheet can be specified:

```
input:
  ont:
    sheet: "test_data/fullPipeline/ont.tsv"
```

The sample sheet has the following content must have the columns `SAMPLE` and `READS` where `SAMPLE` is the sample name and `READS` points to the input reads.
`READS` can be local paths, URLs or S3 paths. In case S3 is used, additional [configuration](pipeline_input.md#s3-Configuration) is necessary.

Example:
```
---8<--- "test_data/fullPipeline/ont.tsv"
```

#### Command line

For a few number of samples it is also possible to specify them on the command line:  

```BASH
	  --input.ont.r read1.fq.gz  \
	  --input.ont.names test1
```

where

  * The `--input.ont.r` parameter points to the Nanopore reads. It uses the same type of resources (URL, S3, etc.) as the `READS` column (see the "Sample Sheet" section).

  * `--input.ont.names` are the sample names. If multiple samples are provided they must be enclosed in double ticks (e.g. "test1 test2").
    The number of sample names must match the number of files provided for `--input.ont.r`.


### Sequence Read Archive (SRA) Samples

The toolkit is able to fetch fastq files based on SRA run accession IDs from the NCBI or from a mirror based on S3:

#### SRA Mirror

!!! info "Additional Configuration"
    For accessing the mirror via S3, additional [configuration](pipeline_input.md#s3-Configuration) is necessary.

For the SRA Mirror it is possible to specify a SRA run accession or SRA project ID.

##### Sample Sheet

```YAML
input:
  SRA:
    pattern:
      ont: ".+[^(_1|_2)].+fastq.gz$"
      illumina: ".+(_1|_2).+fastq.gz$"
    S3:
      sheet: test_data/SRA/ONTsamples.tsv 
      bucket: "s3://ftp.era.ebi.ac.uk" 
      prefix: "/vol1/fastq/"
      watch: false
      skipDB: false
```

where:
  * `sheet` is the path to a file containing a column with `ACCESSION` as header. The `ACCESSION` column contains either SRA run or study accessions.

  * `bucket` is the S3 Bucket hosting the data.

  * `prefix` is the path to the actual SRA datasets.

  * `watch` if true, the file specified with the `sheet` attribute is watched and every time a new SRA run ID is
     appended, the pipeline is triggered. The pipeline will never finish in this mode. Please note that watch currently only works
     if only one input type is specified (e.g "ont" or "paired" ...)
 
  *  `patternONT` and `patternIllumina` are patterns that are applied on the specified mirror in order to select the correct input files.

  *  `skipDB` is a flag for skipping SRA run IDs from the local SRA database. This flag is only used for debugging.

The sample sample sheet must have the column name `ACCESSION`.

Example:
```
---8<--- "test_data/SRA/ONTsamples.tsv"
```

##### Command line

Rather than creating a sample sheet for processing a few samples, you can also specify the SRA run accessions on the command line.

Example:
```BASH
	  --input.SRA.S3.id "SRR29912082 ERR12263778"  
```

#### NCBI SRA

With the following mode SRA datasets can directly be fetched from NCBI.

##### Sample Sheet

```YAML
input:
  SRA:
    pattern:
      ont: ".+[^(_1|_2)].+fastq.gz$"
      illumina: ".+(_1|_2).+fastq.gz$"
    NCBI:
      sheet: test_data/SRA/ncbi_samples.tsv 
```

The sample sample sheet must have the column name `ACCESSION`.

Example:
```
---8<--- "test_data/SRA/ONTsamples.tsv"
```

##### Command line

You can also specify the SRA run accessions on the command line.

Example:
```BASH
	  --input.SRA.NCBI.id "SRR29912082 ERR12263778"  
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
