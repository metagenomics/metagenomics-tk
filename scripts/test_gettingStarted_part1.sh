# This file will be referenced on the online wiki
NXF_HOME=$PWD/.nextflow NXF_VER=23.10.0 nextflow run metagenomics/metagenomics-tk -work-dir $(pwd)/work \
    -profile slurm \
    -entry wFullPipeline \
    -params-file  https://raw.githubusercontent.com/metagenomics/metagenomics-tk/refs/heads/feat/getting_started_guide/default/fullPipeline_illumina_nanpore_getting_started_part1.yml \
    --s3SignIn false \
    --scratch /vol/scratch \
    --databases /vol/scratch/databases \
    --output output \
    --input.paired.path https://raw.githubusercontent.com/metagenomics/metagenomics-tk/refs/heads/feat/getting_started_guide/test_data/fullPipeline/reads_split.tsv

