# Assembly

## Input

=== "Command for short read data"

    ```
    -entry wShortReadAssembly -params-file example_params/assembly.yml
    ```

=== "Command for long read data"

    ```
    -entry wOntAssembly -params-file example_params/assemblyONT.yml
    ```

=== "Megahit Configuration File"

    ```YAML
    ---8<--- "../example_params/assembly.yml"
    ```

=== "Metaspades Configuration File"

    ```YAML
    ---8<--- "../example_params/assemblyMetaspades.yml"
    ```

=== "MetaFlye Configuration File"

    ```YAML
    ---8<--- "../example_params/assemblyONT.yml"
    ```

=== "TSV Table Short Read"

    ```TSV
    ---8<--- "../test_data/assembly/samples.tsv"
    ```

=== "TSV Table Nanopore"

    ```TSV
    ---8<--- "../test_data/assembly/samplesONT.tsv"
    ```
 
 
## Output

The output is a gzipped fasta file containing contigs.

## Megahit Error Handling

On error with exit codes ([-9, 137]) (e.g. due to memory restrictions), the tool is executed again with higher cpu and memory values.
The memory and cpu values are computed by the formula 2^(number of attempts) * (cpu/memory value of the assigned flavour).
The highest possible cpu/memory value is restricted by the highest cpu/memory value of all flavours defined in the resource section 
(see global [configuration](../pipeline_configuration.md) section). 

