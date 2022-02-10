# Annotation

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

Diamond and KeGGFromDiamond need a specific diamond and kegg database as input.
See [database section](##-Database-input-configuration) for possible download strategies.
The compressed kegg database must be a tar.gz file and the diamond file must gzipped. 

## Output

### Diamond

Calculated significant matches of a nucleotide/protein query which was compared against a database.

### Prodigal

Predicted genes in the `*.faa` `*.fna` and `*.gff` file format.

### KEGGFromDiamond

Result `*.tsv` file filled with KEGG informations (linke modules, KO's, ...) which could be linked to the input Diamond hits.
  
### Resistance Gene Identifier (rgi)

The `*rgi.tsv` files contain the found CARD genes.


