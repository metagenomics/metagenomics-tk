# MagAttributes

## Input


=== "Command"

    ```
    -entry wMagAttributes -params-file example_params/magAttributes.yml 
    ```

=== "Configuration File"

    !!! warning "Warning"
     
        **The configuration file shown here is for demonstration and testing purposes only. 
          Parameters that should be used in production can be viewed in the magAttributes section 
          of one of the yaml files located in the `default` folder of the Toolkit's Github repository.**

    ```YAML
    ---8<--- "example_params/magAttributes.yml"
    ```

=== "MAGs TSV Table"

    ```TSV
    ---8<--- "test_data/magAttributes/input.tsv"
    ```
    Must include at least `DATASET` identifier and mag specific `PATH` and `BIN_ID` column.

## Databases

Checkm and GTDB need their databases as input. See [database section](../database.md/#database-input-configuration) for possibly download strategies.
The GTDB and Checkm compressed databases must be tar.gz files. If you provide the extracted version of GTDB using the `extractedDBPath` parameter,
please specify the path to the `releasesXXX` directory (e.g. "/vol/spool/gtdb/release202").

If you need credentials to access your files via S3 then please use the following command:

For GTDB:
```
nextflow secrets set S3_gtdb_ACCESS XXXXXXX
nextflow secrets set S3_gtdb_SECRET XXXXXXX
```

For Checkm:
```
nextflow secrets set S3_checkm_ACCESS XXXXXXX
nextflow secrets set S3_checkm_SECRET XXXXXXX
```

## Output

### GTDBTk

All GTDB files include the GTDB specific columns in addition to a `SAMPLE` column (`SAMPLE_gtdbtk.bac120.summary.tsv`, `SAMPLE_gtdbtk.ar122.summary.tsv`).
In addition, this module produces a file `SAMPLE_gtdbtk_CHUNK.tsv` that combines both files and adds a `BIN_ID` column that adheres to the magAttributes specification

### Checkm and Checkm2

The Checkm and Checkm2 output adheres to the magAttributes specification and adds a `BIN_ID` and `SAMPLE` column to the output file.
If Checkm2 and Checkm are both specified in the config file then only the Checkm2 results are used for downstream pipeline steps.

