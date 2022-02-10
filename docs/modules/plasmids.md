# Plasmids

The plasmid module is able to identify contigs as plasmids and it is able to assemble plasmids from the samples fastq data.
For running the plasmid assembly we suggest to run the full pipeline mode with the enabled plasmids module. See input example configuration files section.

## Input

=== "Command"

    ```
    -entry wPlasmidsPath -params-file example_params/plasmids.yml
    ```

=== "Configuration file for full pipeline mode with plasmids detections"

    ```YAML
    ---8<--- "../example_params/fullPipeline_fraction/fullPipeline_fraction_plasmid.yml"
    ```

=== "Configuration file for plasmids module only"

    ```YAML
    ---8<--- "../example_params/plasmid.yml"
    ```

=== "TSV Table"

    ```TSV
    ---8<--- "../test_data/plasmid/input_contigs.tsv"
    ```

### Databases

PLSDB needs a plasmid database as input. See [database section](##-Database-input-configuration) for possible download strategies.
The compressed database must be a zip file. 

## Output

### SCAPP

SCAPP detects plasmid sequences out of the samples assembly graph.
It reports sequences as gzipped fasta files (`*_plasmids.fasta.gz`). A basic statistic (`*_plasmids_stats.tsv`)
is also generated.

### PlasClass

PlasClass is able to identify contigs as plasmids. It reports gzipped fata files (`*_plasmids.fasta.gz`)
and their probabilities (`*_probs.tsv`).

### PLSDB

PLSDB includes a curated set of plasmid sequences that were extracted from databases like refseq.
The metadata of found sequences are reported in `*.tsv` and the metadata of the filtered sequences in `*_kmerThreshold_X.tsv`.

PLSDB needs a plasmid database as input. See [database section](##-Database-input-configuration) for possibly download strategies.
The compressed database must be a zip file. 

