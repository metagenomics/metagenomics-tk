# Run the following command to test the MetaSpades assembler 
NXF_HOME=$PWD/.nextflow NXF_VER=23.10.0 nextflow run metagenomics/metagenomics-tk \ 
    -work-dir /vol/spool/work \
    -profile slurm \
    -entry wAggregatePipeline \
    -params-file  https://raw.githubusercontent.com/metagenomics/metagenomics-tk/refs/heads/feat/getting_started_guide/default/fullPipelineAggregate.yml \
    --s3SignIn false \
    --scratch /vol/scratch \
    --databases /vol/scratch/databases \
    --input output \
    --output output
