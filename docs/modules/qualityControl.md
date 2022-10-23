# Quality Control

The quality control module removes adapters, trims and filters short read and long read data.

## Input

=== "Command for short read data"

    ```
    -entry wShortReadQualityControl -params-file example_params/qc.yml
    ```

=== "Command for nanopore data"

    ```
    -entry wOntQualityControl -params-file example_params/qcONT.yml
    ```


=== "TSV Table short read"

    ```TSV
    ---8<--- "../test_data/qc/reads_split.tsv"
    ```

=== "TSV Table nanopore"

    ```TSV
    ---8<--- "../test_data/qcONT/ont.tsv"
    ```
 
## Output

The output is a gzipped fastq file (short read: `SAMPLE_interleaved.qc.fq.gz`, long read: `SAMPLE_qc.fq.gz`)
containing trimmed and quality filtered reads.
