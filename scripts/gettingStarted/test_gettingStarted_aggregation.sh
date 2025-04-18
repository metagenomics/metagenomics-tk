# Run the following command to test the MetaSpades assembler 
NXF_HOME=$PWD/.nextflow NXF_VER=24.10.4 nextflow run metagenomics/metagenomics-tk \
    -work-dir $(pwd)/work \
    -profile slurm \
    -ansi-log false \
    -entry wAggregatePipeline \
    -params-file  https://raw.githubusercontent.com/metagenomics/metagenomics-tk/refs/heads/master/default/fullPipelineAggregate.yml \
    --logDir logAggregate \
    --s3SignIn false \
    --scratch /vol/scratch \
    --databases /vol/scratch/databases \
    --input my_data_spades_output \
    --output output

# Check the formed species cluster
cat  my_data_spades_output/AGGREGATED/1/dereplication/*/bottomUpClustering/clusters/clusters.tsv


make check LOG_DIR=$(pwd)/logAggregate
