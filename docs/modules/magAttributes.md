# MagAttributes

## Input


=== "Command"

    ```
    -entry wMagAttributes -params-file example_params/magAttributes.yml 
    ```

=== "Configuration File"

    ```YAML
    ---8<--- "../example_params/magAttributes.yml"
    ```

=== "MAGs TSV Table"

    ```TSV
    ---8<--- "../test_data/magAttributes/input.tsv"
    ```
    Must include at least `DATASET` identifier and mag specific `PATH` and `BIN_ID` column.

## Output

### Prokka

Prokka computes `*.err`, `*.faa`, `*.ffn`, `*.fna`, `*.fsa`, `*.gbk`, `*.gff`, `*.sqn`, `*.tbl`, `*.tbl` for every bin.
Details of all files can be read on the Prokka page.
In addition it also computes a summary tsv file which adheres to the magAttributes specification.

### GTDBTk

All GTDB files include the GTDB specific columns in addition to a `SAMPLE` column (`SAMPLE_gtdbtk.bac120.summary.tsv`, `SAMPLE_gtdbtk.ar122.summary.tsv`).
In addition, this module produces a file `SAMPLE_gtdbtk_CHUNK.tsv` that combines both files and adds a `BIN_ID` column that adheres to the magAttributes specification

### Checkm

The Checkm output adheres to the magAttributes specification and adds a `BIN_ID` and `SAMPLE` column to the output file.


