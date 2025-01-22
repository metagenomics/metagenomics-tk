# Export

This module exports a set of results produced by Metagenomics-tk.
Currently, the export for EMGB needs the results of the following tools:

1. Assembly
2. Binning
3. Checkm(2)
4. Prokka output
5. gtdb-tk 
6. Mmseqs Taxonomy (Database: GTDB)
7. MMseqs (Database: uniref90) 

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

* blastDB: The toolkit runs mmseqs against multiple databases. You can specify here, which blast output should be used. (Default: uniref90)

* taxonomyDB: MMseqs is executed against a specific taxonomy database. (Default: gtdb)   

## Output

The following files are produced as output:

* SAMPLE.bins.json.gz  
* SAMPLE.contigs.json.gz
* SAMPLE.genes.json.gz

where `SAMPLE` is the name of the sample.

You can read [here](https://gitlab.ub.uni-bielefeld.de/cmg/emgb/emgb-server/-/tree/master#quick-start) how to 
start EMGB and how to use these files to import a dataset. 
