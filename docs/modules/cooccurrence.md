# Cooccurrence

```
-entry wCooccurrence -params-file example_params/coocurrence.yml
```

## Input

=== "Command"

    ```
    -entry wCooccurrence -params-file example_params/coocurrence.yml
    ```

=== "Configuration File"

    ```YAML
    ---8<--- "../example_params/coocurrence.yml"
    ```

=== "TSV Table"

    ```TSV
    ---8<--- "../test_data/assembly/samples.tsv"
    ```
    Contains abundance values of mags per sample.
 
=== "GTDB TSV Table"

    ```TSV
    ---8<--- "../test_data/cooccurrence/gtdb.tsv"
    ```
    GTDB assignmend of all samples that were produced by magAttributes module.

 
## Output

 * Graphml file for further processing 


