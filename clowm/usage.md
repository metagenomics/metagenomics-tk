# Meta-Omics-Toolkit: Usage

### Paired End Input
The input should be a path to a tsv file containing a sample id, as well as a path to the left and right read.

```
SAMPLE  READS1  READS2
test1 s3://my-bucket/path/to/read1_1.fq.gz  s3://my-bucket/path/to/read1_2.fq.gz
test2 s3://my-bucket/path/to/read2_1.fq.gz  s3://my-bucket/path/to/read2_2.fq.gz
```

### Nanopore Input
For Nanopore data a seperate input file should be specified.

```
SAMPLE  READS
nano  s3://my-bucket/path/to/reads.fq.gz
```
