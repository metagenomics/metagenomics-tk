# Run the following command to test the MetaSpades assembler 
NXF_HOME=$PWD/.nextflow NXF_VER=23.10.0 nextflow run metagenomics/metagenomics-tk \ 
    -work-dir /vol/spool/work \
    -profile slurm \
    -ansi-log false \
    -entry wAggregatePipeline \
    -params-file  https://raw.githubusercontent.com/metagenomics/metagenomics-tk/refs/heads/master/default/fullPipelineAggregate.yml \
    --s3SignIn false \
    --scratch /vol/scratch \
    --databases /vol/scratch/databases \
    --input my_data_spades_output \
    --output output

# Check the formed species cluster
cat  my_data_spades_output/AGGREGATED/1/dereplication/*/bottomUpClustering/clusters/clusters.tsv
