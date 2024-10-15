# Quality Control

The quality control module removes adapters, trims and filters short read and long read data.

## Input

=== "Command for short read data"

    ```BASH
    -entry wShortReadQualityControl -params-file example_params/qc.yml
    ```

=== "Configuration File"

    !!! warning "Warning"
     
        **The configuration file shown here is for demonstration and testing purposes only. 
          Parameters that should be used in production can be viewed in the quality control section 
          of one of the yaml files located in the `default` folder of the Toolkit's Github repository.**

    ```YAML
    ---8<--- "example_params/qc.yml"
    ```

=== "Command for nanopore data"

    ```
    -entry wOntQualityControl -params-file example_params/qcONT.yml
    ```


=== "TSV Table short read"

    ```TSV
    ---8<--- "test_data/qc/reads_split.tsv"
    ```

=== "TSV Table nanopore"

    ```TSV
    ---8<--- "test_data/qcONT/ont.tsv"
    ```
 
## Output

The output is a gzipped fastq file (short read: `SAMPLE_interleaved.qc.fq.gz`, long read: `SAMPLE_qc.fq.gz`)
containing trimmed and quality filtered reads.
