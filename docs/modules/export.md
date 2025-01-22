# Export

This module exports a set of results produced by Metagenomics-tk.
Currently, the export for EMGB needs the results of the following tools:

1. Assembly
2. Binning
3. Checkm (v1 or v2)
4. Prokka output
5. GTDB-Tk 
6. MMseqs Taxonomy (Database: GTDB)
7. MMseqs (Database: UniRef90) 

## Input

=== "Command"

    ```
    -entry wExportPipeline -params-file example_params/export.yml
    ```

=== "Configuration File"

    !!! warning "Warning"
     
        **The configuration file shown here is for demonstration and testing purposes only. 
          Parameters that should be used in production can be viewed in the read mapping section 
          of one of the yaml files located in the `default` folder of the Toolkit's Github repository.**

    ```YAML
    ---8<--- "example_params/export.yml"
    ```


### Additional Parameters

* blastDB: The toolkit runs MMseqs against multiple databases. You can specify here, which BLAST output should be used. (Default: UniRef90)

* taxonomyDB: MMseqs is executed against a specific taxonomy database. (Default: GTDB)   

## Output

The following files are produced as output:

* SAMPLE.bins.json.gz  
* SAMPLE.contigs.json.gz
* SAMPLE.genes.json.gz

where `SAMPLE` is the name of the sample.

You can read [here](https://gitlab.ub.uni-bielefeld.de/cmg/emgb/emgb-server/-/tree/master#quick-start) how to 
start EMGB and how to use these files to import a dataset. 
