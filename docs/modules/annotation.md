# Annotation

The annotation module is able to predict genes and annotate those based on a set of user provided databases.
A user can add additional Diamond formatted databases as part of the diamond configuration by adding a key (Example: `kegg` ) with 
a possible download strategy. See [database section](##-Database-input-configuration) for possible download strategies.
In addition, the resistance gene identifier is executed by default.

## Input  

=== "Command"

    ```
    -entry wAnnotateLocal -params-file example_params/annotation.yml
    ```

=== "Configuration File"

    ```YAML
    ---8<--- "../example_params/annotation.yml"
    ```

=== "TSV Table"

    ```TSV
    ---8<--- "../test_data/annotation/input_small.tsv"
    ```

### Databases

KeGGFromDiamond is only executed if genes are searched against the kegg diamond database. In the diamond section there must be a diamond run specified with `kegg` as 
identifier (see example configuration file).
KeGGFromDiamond need a kegg database as input which must be a tar.gz file.
See [database section](##-Database-input-configuration) for possible download strategies.

## Output

### Diamond

Calculated significant matches of a nucleotide/protein query which was compared against a user provided set of databases.

### Prodigal

Predicted genes in the `*.faa` `*.fna` and `*.gff` file format.

### KEGGFromDiamond

Result `*.tsv` file filled with KEGG informations (linke modules, KO's, ...) which could be linked to the input Diamond hits.
  
### Resistance Gene Identifier (rgi)

The `*rgi.tsv` files contain the found CARD genes.


