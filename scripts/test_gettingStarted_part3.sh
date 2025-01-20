# Run the following command with files stored on local disk
NXF_HOME=$PWD/.nextflow NXF_VER=23.10.0 nextflow run main.nf -work-dir $(pwd)/work \
    -profile slurm \
    -entry wFullPipeline \
    -params-file https://github.com/metagenomics/metagenomics-tk/blob/master/default/fullPipeline_illumina_nanpore_getting_started_part3.yml
    --s3SignIn false \
    --scratch /vol/scratch \
    --databases /vol/scratch/databases \
    --output _data_output \
    --input.paired.path inputFiles/input.tsv

