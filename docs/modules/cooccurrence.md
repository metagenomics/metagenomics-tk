# Cooccurrence

The Cooccurrence module builds a cooccurrence network where each node is a MAG
and every edge represents an association between them. The network can be inferred based on
correlation or inverse covariance estimation by SPIEC-EASI. SPIEC-EASI is executed multiple times
based on different parameter settings in order to find the most stable network.
In addition, it is possible to compute multiple metrics for every edge based on genome-scale
metabolic models and the SMETANA metrics. 

```
-entry wCooccurrence -params-file example_params/cooccurrence.yml
```

## Input

=== "Command"

    ```
    -entry wCooccurrence -params-file example_params/cooccurrence.yml
    ```
=== "Configuration File for Cooccurrence"

    !!! warning "Warning"
     
        **The configuration file shown here is for demonstration and testing purposes only. 
          Parameters that should be used in production can be viewed in the cooccurrence section 
          of one of the yaml files located in the `default` folder of the Toolkit's Github repository.**

    ```YAML
    ---8<--- "example_params/cooccurrence.yml"
    ```

=== "TSV Table"

    ```TSV
    ---8<--- "test_data/assembly/samples.tsv"
    ```
    Contains abundance values of mags per sample.
 
=== "GTDB TSV Table"

    ```TSV
    ---8<--- "test_data/cooccurrence/gtdb.tsv"
    ```
    GTDB assignment of all samples that were produced by magAttributes module.

=== "Configuration File for analyzing edges in Cooccurrence Graph"

    !!! warning "Warning"
     
        **The configuration file shown here is for demonstration and testing purposes only. 
          Parameters that should be used in production can be viewed in the cooccurrence section 
          of one of the yaml files located in the `default` folder of the Toolkit's Github repository.**

    ```YAML
    ---8<--- "example_params/cooccurrenceMetabolom.yml"
    ```

=== "GTDB TSV for analyzing Edges"

    ```TSV
    ---8<--- "test_data/cooccurrence/gtdb_large.tsv"
    ```
    GTDB assignment of all samples that were produced by the magAttributes module.


=== "Model TSV for computing Metabolomics Metrics on Edges"

    ```TSV
    ---8<--- "test_data/cooccurrence/models.tsv"
    ```

The following parameters can be configured:
  
  * metabolicEdgeBatches: Batches of edges that are provided as input to SMETANA.
     
  * metabolicEdgeReplicates: Number of replicates per edge that should be computed.

## Output

 * output_raw.graphml: Cooccurrence network unfiltered in graphml format

 * output.graphml: Filtered cooccurrence network in graphml format

 * edges_index.tsv: Edges of the graph

 * edgeAttributes.tsv: Edge attributes of the graph containing metrics computed by SMETANA.

### SPIEC-EASI

 * stability.txt: Network stability estimation

