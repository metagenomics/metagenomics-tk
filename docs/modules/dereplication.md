# Dereplication

## Input

=== "Command"

    ```
    -entry wDereplication -params-file example_params/dereplication.yml
    ```

=== "Configuration File"

    ```YAML
    ---8<--- "../example_params/dereplication.yml"
    ```

=== "TSV Table"

    ```TSV
    ---8<--- "../test_data/dereplication/input.tsv"
    ```
    Must include the columns `DATASET`, `BIN_ID`, `PATH`, `COMPLETENESS`, `CONTAMINATION`, `COVERAGE`, `N50` and `HETEROGENEITY`. 
    Completeness and contamination can be used for filtering (see `params-file`). `N50`, `COVERAGE` and `HETEROGENEITY` are used for selecting the representative of every cluster.
    You can set values of these columns to zero if data is not available or if you don't want the representative selection to be influenced by theses columns. Make sure that `BIN_ID` is a unique identifier.

## Output

The output tsv file (`final_clusters.tsv`) contains the columns `CLUSTER`, `GENOME` and `REPRESENTATIVE` where `CLUSTER` identifies a group of genomes, `GENOME` represents the path or
link of a genome and `REPRESENTATIVE` is either 0 or 1 (selected as representative).
If `sans` is specified in the configuration file (see examples folder), then [SANS](https://gitlab.ub.uni-bielefeld.de/gi/sans) is used to dereplicate every cluster reported by the previous step further down 
to generate strain specific clusters. 

