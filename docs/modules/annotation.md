# Annotation

The annotation module is able to predict genes and annotate those based on a set of user provided databases.
A user can add additional formatted databases as part of the configuration by adding a key (Example: `kegg` ) with 
a possible download strategy. See [database section](../pipeline_configuration.md##-Database-input-configuration) for possible download strategies.
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

KeGGFromBlast is only executed if genes are searched against a KEGG database. There must be a `kegg` identifier (see example configuration file) in the annotation section.
KeGGFromBlast needs a kegg database as input which must be a tar.gz file.
See [database section](../pipeline_configuration.md##-Database-input-configuration) for possible download strategies.

## Output

### Diamond

Calculated significant matches of a nucleotide/protein query which was compared against a user provided set of databases.

### MMseqs2

Calculated significant matches of a nucleotide/protein query which was compared against a user provided set of databases.

### MMseqs2 - Taxonomy

By identifying homologs through searches against a provided MMseqs2 taxonomy-database, MMseqs2 can compute the lowest common ancestor. 
This lowest common ancestor is a robust taxonomic label for unknown sequences.
These labels are presented in form of an `*.taxonomy.tsv` file, a `*.krakenStyleTaxonomy.out` formated in accordance to the KRAKEN tool outputs and
an interactive KRONA plot in form of a html website `*.krona.html`.

### Prodigal

Predicted genes in the `*.faa` `*.fna` and `*.gff` file format.

### KEGGFromBlast

Result `*.tsv` file filled with KEGG information (like modules, KO's, ...) which could be linked to the input hits.
  
### Resistance Gene Identifier (rgi)

The `*rgi.tsv` files contain the found CARD genes.


