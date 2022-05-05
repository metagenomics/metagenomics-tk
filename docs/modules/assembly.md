# Assembly

## Input

=== "Command"

    ```
    -entry wAssembly -params-file example_params/assembly.yml
    ```

=== "Megahit Configuration File"

    ```YAML
    ---8<--- "../example_params/assembly.yml"
    ```

=== "Metaspades Configuration File"

    ```YAML
    ---8<--- "../example_params/assemblyMetaspades.yml"
    ```

=== "TSV Table"

    ```TSV
    ---8<--- "../test_data/assembly/samples.tsv"
    ```
 
## Output

The output is a gzipped fasta file containing contigs.

## Megahit Error Handling

On error with exit codes ([-9, 137]) (e.g. due to memory restrictions), the tool is executed again with higher cpu and memory values.
The memory and cpu values are computed by the formula 2^(number of attempts) * (cpu/memory value of the assigned flavour).
The highest possible cpu/memory value is restricted by the highest cpu/memory value of all flavours defined in the resource section 
(see global [configuration](../pipeline_configuration.md) section). 

