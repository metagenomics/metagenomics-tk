
# Module Specification

## Assembly

* Version: 0.2.0

### Output:

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
## Binning

* Version: 0.2.0

### Output:

Binning file names must fulfill the following name pattern:

```
SAMPLENAME_bin.NUMBER.fa
```

Where `NUMBER` is a unique identifier per SAMPLE.

Contig names must be renamed according to the following pattern:

`SAMPLEID_SEQUENCECOUNTER_SEQUENCEHASH MAG=BINNUMBER`

where

   * `SAMPLEID`, `SEQUENCECOUNTER` and `SEQUENCEHASH`  definitions can be inspected in the assembly specification.

   * `BINNUMBER` is an unique identifier per SAMPLE.

## MAG Attributes

* Version: 0.1.0

### Output:

```
SAMPLENAME_TOOLNAME_CHUNK.tsv
```

where
 * `TOOLNAME` could be for example `checkm`, `gtdb` etc.
 * `CHUNK` is a random identifier that produces values for one part of all MAGs of a given sample.


### FORMAT

The header line specifies the following columns: 

```
BIN_ID	SAMPLE	BIN_ATTRIBUTE1	BIN_ATTRIBUTE2 ...
```

where 
  * `BIN_ID` is unique for the samples.
  * `BIN_ATTRIBUTES` All column names that are not `BIN_ID` or `SAMPLE` can be any property of a MAG, like contamination, completeness etc.

## Quality Control

* Version: 0.1.0

### Output:

```
SAMPLE_interleaved.fq.gz
```

