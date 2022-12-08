# Annotation

The annotation module is able to predict genes and annotate those based on Prokka and a set of user provided databases.
A user can add additional formatted databases as part of the configuration by adding a key (Example: `kegg` ) with 
a possible download strategy. See [database section](../pipeline_configuration.md#database-input-configuration) for possible download strategies.
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

### MMseqs2

MMseqs2 needs a combination of different data, index and dbtype files as "one" database, be it in- or output.
See [MMseqs2 database](https://github.com/soedinglab/mmseqs2/wiki#mmseqs2-database-format) for more information.
As multiple and in most cases, big files are used, tar and [zstd](https://github.com/facebook/zstd) are utilized to compress and transport files.
Input databases have to be compressed by these and need to end with `.tar.zst`. Naming inside an archive is irrelevant, as databases are picked automatically.
Multiple databases per one archive are not supported, one archive, one database.

#### KEGGFromBlast
KeGGFromBlast is only executed if genes are searched against a KEGG database. There must be a `kegg` identifier (see example configuration file) in the annotation section.
KeGGFromBlast needs a kegg database as input which must be a tar.gz file.
See [database section](../pipeline_configuration.md#database-input-configuration) for possible download strategies.

### RGI

RGI needs a CARD database which can be fetched via this link:  https://card.mcmaster.ca/latest/data.
The compressed database must be a tar.bz2 file. 
See [database section](../pipeline_configuration.md#database-input-configuration) for possible download strategies.

## Output

### MMseqs2

Calculated significant matches of a nucleotide/protein query which was compared against a user provided set of databases.

### MMseqs2 - Taxonomy

By identifying homologous through searches against a provided MMseqs2 taxonomy-database, MMseqs2 can compute the lowest common ancestor. 
This lowest common ancestor is a robust taxonomic label for unknown sequences.
These labels are presented in form of an `*.taxonomy.tsv` file, a `*.krakenStyleTaxonomy.out` formatted in accordance to the [KRAKEN tool](https://ccb.jhu.edu/software/kraken/) outputs and
an interactive [KRONA](https://github.com/marbl/Krona/wiki) plot in form of a html website `*.krona.html`.

### Prokka

Prokka computes `*.err`, `*.faa`, `*.ffn`, `*.fna`, `*.fsa`, `*.gbk`, `*.gff`, `*.sqn`, `*.tbl`, `*.tbl` for every bin.
Details of all files can be read on the Prokka page.
In addition, it also computes a summary tsv file which adheres to the magAttributes specification.

### KEGGFromBlast

Result `*.tsv` file filled with KEGG information (like modules, KO's, ...) which could be linked to the input hits.
  
### Resistance Gene Identifier (rgi)

The `*rgi.tsv` files contain the found CARD genes.


