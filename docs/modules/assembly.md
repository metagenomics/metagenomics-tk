# Assembly

## Input

=== "Command for short read data with optional single end reads"

    ```
    -entry wShortReadAssembly -params-file example_params/assembly.yml
    ```

=== "Command for long read data"

    ```
    -entry wOntAssembly -params-file example_params/assemblyONT.yml
    ```

=== "Megahit Configuration File"

    !!! warning "Warning"
     
        **The configuration file shown here is for demonstration and testing purposes only. 
          Parameters that should be used in production can be viewed in the assembly section 
          of one of the yaml files located in the `default` folder of the Toolkit's Github repository.**

    ```YAML
    ---8<--- "example_params/assembly.yml"
    ```

=== "Metaspades Configuration File"

    ```YAML
    ---8<--- "example_params/assemblyMetaspades.yml"
    ```

=== "MetaFlye Configuration File"


    !!! warning "Warning"
     
        **The configuration file shown here is for demonstration and testing purposes only. 
          Parameters that should be used in production can be viewed in the assemblyONT section 
          of one of the yaml files located in the `default` folder of the Toolkit's Github repository.**

    ```YAML
    ---8<--- "example_params/assemblyONT.yml"
    ```

=== "TSV Table Short Read"

    ```TSV
    ---8<--- "test_data/assembly/samples.tsv"
    ```

=== "TSV Table Nanopore"

    ```TSV
    ---8<--- "test_data/assembly/samplesONT.tsv"
    ```
 
 
## Output

The output is a gzipped fasta file containing contigs.

## Megahit 

### Error Handling

On error with exit codes ([-9, 137, 247]) (e.g. due to memory restrictions), the tool is executed again with higher cpu and memory values.
The memory and cpu values are in case of a retry selected based on the flavor with the next higher memory value.
The highest possible cpu/memory value is restricted by the highest cpu/memory value of all flavors defined in the resource section 
(see global [configuration](../configuration.md) section). 

### Peak memory usage prediction

Memory consumption of an assembler varies based on diversity. We trained a machine learning model on kmer frequencies
and the nonpareil diversity index in order to be able to predict the memory peak consumption of megahit in our full pipeline mode. The required
resources in order to run the assembler are thereby fitted to the resources that are actually needed for a specific dataset. If this
mode is enabled then Nonpareil and kmc that are part of the quality control module are automatically executed before the assembler run.  

Please note that this mode is only tested for Megahit with default parameters and the meta-sensitive mode (`--presets meta-sensitive`).

```
  resources:
    RAM: 
      mode: MODE
      predictMinLabel: LABEL
```

where 
    * MODE can be either 'PREDICT' for predicting memory usage or 'DEFAULT' for using a default flavor defined in the resources section.

    * LABEL is the flavor that will be used if the predicted RAM is below the memory value defined as part of the LABEL flavor. It can also be set to AUTO to always use the predicted flavor.
