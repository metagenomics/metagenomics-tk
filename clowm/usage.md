# Metagenomics-Toolkit: Usage

### Paired End Input
The input should be a path to a tsv file containing a sample id, as well as a link to the left and right read. If your tools are publicly hosted 
on a server, then you can specify an HTTPs link. If you have uploaded fastq files to a bucket in the CloWM service then you can specify the 
bucket name as part of the s3 link (see example).

#### Paired End Input Bucket Example
```
SAMPLE  READS1  READS2
test1 s3://my-bucket/path/to/read1_1.fq.gz  s3://my-bucket/path/to/read1_2.fq.gz
test2 s3://my-bucket/path/to/read2_1.fq.gz  s3://my-bucket/path/to/read2_2.fq.gz
```

where 
  * "my-bucket" is the bucket created in the CloWM service.
  * "/path/to" is the the path to the fastq files in the CloWM bucket.

#### Paired End Input HTTPS Example

```
SAMPLE  READS1  READS2
test1   https://openstack.cebitec.uni-bielefeld.de:8080/swift/v1/meta_test/small/read1_1.fq.gz  https://openstack.cebitec.uni-bielefeld.de:8080/swift/v1/meta_test/small/read2_1.fq.gz
test2   https://openstack.cebitec.uni-bielefeld.de:8080/swift/v1/meta_test/small/read1_1.fq.gz  https://openstack.cebitec.uni-bielefeld.de:8080/swift/v1/meta_test/small/read2_1.fq.gz
```

All fastq files in this example are publicly available.