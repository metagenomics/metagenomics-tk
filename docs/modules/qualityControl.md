# Quality Control

The quality control module removes adapters, trims and filters short and read long read data.
Since quality control is the usual first step in the processing of sequencing data, the Toolkit offers a way to
directly download the sequencing data (See `download` flag.). This allows the data to be downloaded in parallel on multiple machines,
as opposed to the usual Nextflow mechanism of downloading input data only on the VM running Nextflow.
In addition, the quality control module allows human reads to be filtered out and, with Nonpareil, provides diversity estimation of input sequences. 

## Short Reads

For short reads we offer a way to only generate a quality report using fastp. This way no additional
disk space for quality controlled reads is necessary. (See `reportOnly` flag in the configuration file below.)

### Input

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

=== "TSV Table short read"

    ```TSV
    ---8<--- "test_data/qc/reads_split.tsv"
    ```

### Output

#### Fastp

`SAMPLE_fastp.json`

:  Contains quality statistics about the raw reads and the quality controlled reads in json format. 

`SAMPLE_fastp_summary_after.tsv`

:  Contains quality statistics about the quality controlled reads in tsv format. 

`SAMPLE_fastp_summary_before.tsv`

:  Contains quality statistics about the raw reads in tsv format. 

`SAMPLE_interleaved.qc.fq.gz`

:  Quality controlled reads

`SAMPLE_report.html`

:  HTML report with plots summarizing the quality of the raw and quality controlled reads.

`test1_unpaired.qc.fq.gz`

: Unpaired reads where the other pair was filtered out due to quality control. 

`test1_unpaired_summary.tsv`

:  TSV file that contains quality statistics about single reads where the other pair was filtered out due to quality control. 

#### KMC 

`SAMPLE.[13|21|71].kmc.json`

K-mer statistics for k-mers of length 13, 21 and 71.

`SAMPLE.[13|21|71].histo.tsv`

K-mer frequency table with the columns `FREQUENCY`, `COUNT` and `SAMPLE`.
`FREQUENCY` is the number of times a specific k-mer appears.
`COUNT` is the number of different k-mers that occur a number of times described by `FREQUENCY`.

## Nanopore Reads

### Input

=== "Command for nanopore data"

    ```
    -entry wOntQualityControl -params-file example_params/qcONT.yml
    ```

=== "Configuration File"

    !!! warning "Warning"
     
        **The configuration file shown here is for demonstration and testing purposes only. 
          Parameters that should be used in production can be viewed in the quality control section 
          of one of the yaml files located in the `default` folder of the Toolkit's Github repository.**

    ```YAML
    ---8<--- "example_params/qcONT.yml"
    ```

=== "TSV Table nanopore"

    ```TSV
    ---8<--- "test_data/qcONT/ont.tsv"
    ```
 
### Output

#### Porechop

`SAMPLE_qc.fq.gz`

: Gzipped quality controlled reads of the format `SAMPLE_qc.fq.gz`.


#### Nanoplot

`NanoStats.tsv`

: Statistics of the output reads, such as quality scores and read length.

`plots`

: Nanoplot offers a variety of plots that show the quality, length and number of the reads.

## Output

The following output is produced for short and long reads. 

### Nonpareil

`SAMPLE.npa`

: You can read more [here](https://nonpareil.readthedocs.io/en/latest/redundancy.html#output)

`SAMPLE.npc`

: You can read more [here](https://nonpareil.readthedocs.io/en/latest/redundancy.html#output)

`SAMPLE.npl`

: You can read more [here](https://nonpareil.readthedocs.io/en/latest/redundancy.html#output)

`SAMPLE_nonpareil_curves.pdf`

: Nonpareil curves visualize the estimated average coverage for the current sequencing effort.

`SAMPLE_nonpareil_index.tsv`

: Nonpareil statistics including the Nonpareil diversity index. 

### Filtered out Human Sequences

`SAMPLE_filtered.fq.gz`

: Sequences without human dna.

`SAMPLE_removed.fq.gz`

: Sequences that were classified as human dna. 

`SAMPLE_summary_[after|before].tsv`

: Statistics of the reads before and after quality control. 
