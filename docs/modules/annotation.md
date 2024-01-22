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
Multiple databases per one archive are not supported, one archive, one database. If the database also includes a taxonomy 
as described [here](https://github.com/soedinglab/mmseqs2/wiki#creating-a-seqtaxdb), it can also be used for taxonomic classifications with MMseqs2 - Taxonomy.
See [database section](../pipeline_configuration.md#database-input-configuration) for possible download strategies.
If you need credentials to access your files via S3 then please use the following command:

```
nextflow secrets set S3_db_ACCESS XXXXXXX
nextflow secrets set S3_db_SECRET XXXXXXX
```

where `db` is the name of the database that you use in your config file.
Example:

```
....
      vfdb:
        params: ' -s 1 --max-seqs 100 --max-accept 50 --alignment-mode 1 --exact-kmer-matching 1 --db-load-mode 3'
        database:
          download:
            source: s3://databases/vfdb_full_2022_07_29.tar.zst
            md5sum: 7e32aaed112d6e056fb8764b637bf49e
            s5cmd:
              params: " --retry-count 30 --endpoint-url https://openstack.cebitec.uni-bielefeld.de:8080 " 
....
```

Based on these settings, you would set the following secret:

```
nextflow secrets set S3_vfdb_ACCESS XXXXXXX
nextflow secrets set S3_vfdb_SECRET XXXXXXX
```

### KEGGFromBlast

KeGGFromBlast is only executed if genes are searched against a KEGG database. There must be a `kegg` identifier (see example configuration file) in the annotation section.
KeGGFromBlast needs a kegg database as input which must be a tar.gz file.
See [database section](../pipeline_configuration.md#database-input-configuration) for possible download strategies.
If you need credentials to access your files via S3 then please use the following command:

```
nextflow secrets set S3_kegg_ACCESS XXXXXXX
nextflow secrets set S3_kegg_SECRET XXXXXXX
```

### MMSeqs Taxonomy

If you need credentials to access your files via S3 then please use the following command:

```
nextflow secrets set S3_TAX_db_ACCESS XXXXXXX
nextflow secrets set S3_TAX_db_SECRET XXXXXXX
```

where `db` is the name of the database that you use in your config file.
Example:

```
....
    mmseqs2_taxonomy:
      gtdb:
        params: ' --orf-filter-s 1 -e 1e-15'
        ramMode: false
        database:
          download:
            source: s3://databases/gtdb_r214_1_mmseqs.tar.gz
            md5sum: 3c8f12c5c2dc55841a14dd30a0a4c718
            s5cmd:
              params: " --retry-count 30 --endpoint-url https://openstack.cebitec.uni-bielefeld.de:8080 " 
....
```

Based on these settings, you would set the following secrets:

```
nextflow secrets set S3_TAX_gtdb_ACCESS XXXXXXX
nextflow secrets set S3_TAX_gtdb_SECRET XXXXXXX
```

### RGI

RGI needs a CARD database which can be fetched via this link:  https://card.mcmaster.ca/latest/data.
The compressed database must be a tar.bz2 file. 
See [database section](../pipeline_configuration.md#database-input-configuration) for possible download strategies.
If you need credentials to access your files via S3 then please use the following command:

```
nextflow secrets set S3_rgi_ACCESS XXXXXXX
nextflow secrets set S3_rgi_SECRET XXXXXXX
```

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
`*.gbk` and `*.sqn` are skipped per default, since tbl2asn runs for quite a while! If you need those files generated by prokka, include:
`--tbl2asn` in the prokka parameters to enable it.
Details of all files can be read on the Prokka page.
In addition, it also computes a summary tsv file which adheres to the magAttributes specification.

### KEGGFromBlast

Result `*.tsv` file filled with KEGG information (like modules, KO's, ...) which could be linked to the input hits.
  
### Resistance Gene Identifier (rgi)

The `*rgi.tsv` files contain the found CARD genes.


