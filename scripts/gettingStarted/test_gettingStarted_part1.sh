# This file will be referenced on the online wiki
NXF_HOME=$PWD/.nextflow NXF_VER=23.10.0 nextflow run metagenomics/metagenomics-tk -work-dir $(pwd)/work \
    -profile slurm \
    -entry wFullPipeline \
    -ansi-log false \
    -params-file  https://raw.githubusercontent.com/metagenomics/metagenomics-tk/refs/heads/master/default/fullPipeline_illumina_nanpore_getting_started_part1.yml \
    --logDir log1 \
    --s3SignIn false \
    --scratch /vol/scratch \
    --databases /vol/scratch/databases \
    --output output \
    --input.paired.path https://raw.githubusercontent.com/metagenomics/metagenomics-tk/refs/heads/master/test_data/fullPipeline/reads_split.tsv


make check LOG_DIR=$(pwd)/logs
