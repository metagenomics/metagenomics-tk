# Assembly

* Version: 0.2.0

## Output:

Assembly file names must fulfill the following name pattern:

```
SAMPLENAME_contigs.fa.gz
```

Contig names must be renamed according to the following pattern:

`SAMPLEID_SEQUENCECOUNTER_SEQUENCEHASH`

where

   * `SAMPLEID` is the name of the dataset (e.g: `SRR234235`)

   * `SEQUENCECOUNTER` is the counter of the contig entry in the fasta file (e.g: 2)

   * `SEQUENCEHASH` are the last 5 characters of an md5sum hash of the fasta entry without the header and newline character. 
     (eg. echo -n "ACGT" | md5sum | cut -d ' ' -f 1 | cut -c -5 )
