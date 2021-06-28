# MAG Attributes

* Version: 0.1.0

## Output:

```
SAMPLENAME_TOOLNAME_CHUNK.tsv
```

where
 * `TOOLNAME` could be for example `checkm`, `gtdb` etc.
 * `CHUNK` is a random identifier that produces values for one part of all MAGs of a given sample.


## FORMAT

The header line specifies the following columns: 

```
BIN_ID	SAMPLE	BIN_ATTRIBUTE1	BIN_ATTRIBUTE2 ...
```

where 
  * `BIN_ID` is of the format `SAMPLE_bin.NUMBER.fa`. NUMBER is unique for the sample.
  * `BIN_ATTRIBUTES` All column names that are not `BIN_ID` or `SAMPLE` can be any property of a MAG, like contamination, completeness etc.
