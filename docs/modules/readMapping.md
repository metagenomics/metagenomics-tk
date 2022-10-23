# Read Mapping

**Note:** This module only supports illumina data. 

## Input

=== "Command"

    ```
    -entry wReadMapping -params-file example_params/readMapping.yml
    ```

=== "Configuration File"

    ```YAML
    ---8<--- "../example_params/readMapping.yml"
    ```

=== "MAGs TSV Table"

    ```TSV
    ---8<--- "../test_data/readMapping/mags.tsv"
    ```

=== "Samples TSV Table"

    ```TSV
    ---8<--- "../test_data/readMapping/samples.tsv"
    ```

## Output

The produced output files are the following: count.tsv, mean.tsv, mean_mincov10.tsv, rpkm.tsv, tpm.tsv, trimmed_mean.tsv.
The content of the files are produced by coverm. All metrics are explained on the coverm GitHub page: https://github.com/wwood/CoverM .


