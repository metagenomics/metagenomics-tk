# Plasmids

The plasmid module is able to identify contigs as plasmids and also to assemble plasmids from the samples fastq data. The module is executed in two
parts. In the first part contigs of a metagenome assembler are scanned for plasmids. In the second part a plasmid assembler is used to assemble
circular plasmids out of raw reads. All plasmid detection tools are executed on the circular assembly result and on the contigs of the metagenome assembler.
But just in the case of the metagenome contigs, the plasmid detection tool results are also used for filtering. Just the filtered contigs are used for downstream
analysis. 

The identification of plasmids is based on the combined result of tools which have a `filter` property assigned. Results of all tools that
have the `filter` property set to true are combined either by a logical `OR` or by a logical `AND`. It is also possible to simply run a tool without
using its result as filter by setting `filter` to `false`. If a tool should be not executed then the tool section should be removed.
Only the detected plasmids will be used for downstream analysis.

All tools that can be used for filtering are also applied for detecting plasmids in circular assemblies.
For running plasmid assembly we suggest to run the full pipeline mode with the enabled plasmids module. See input example configuration files.

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

The plasmid module needs the following compressed database file formats: 

* ViralVerifyPlasmid: pfam.hmm.gz

* MobTyper: mob_suite.tar.gz

* Platon: platon.tar.gz

* PLSDB: plsdb.zip

See [database section](##-Database-input-configuration) for possible download strategies.

#### PLSDB

PLSDB needs a plasmid database as input. See [database section](##-Database-input-configuration) for possible download strategies.
The compressed database must be a zip file. 

## Output

### SCAPP

SCAPP detects plasmid sequences out of the samples assembly graph.
It reports sequences as gzipped fasta files (`*_plasmids.fasta.gz`). A basic statistic (`*_plasmids_stats.tsv`)
is also generated. Coverm coverage metrics are generated for all plasmids. Gene coverage values are generated as part of the annotation module output.

### PlasClass

PlasClass is able to identify plasmids by using a statistical model that was build using kmer frequencies.
It reports gzipped fata files and their probabilities (`*_plasclass.tsv`).

### MobTyper and Platon

MobTyper and Platon are using both replicon typing for plasmid detection. (`*_mobtyper_results.tsv`, `*_platon.tsv`)

### ViralVerifyPlasmid

ViralVerfiy is applying a Naive Bayes classifier (`*_viralverifyplasmid.tsv`).

### PLSDB

PLSDB includes a curated set of plasmid sequences that were extracted from databases like refseq.
The metadata of found sequences are reported in `*.tsv` and the metadata of the filtered sequences in `*_kmerThreshold_X.tsv`.

