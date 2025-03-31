# This file will be referenced on the online wiki

cd ~/mgcourse/

NXF_VER=24.10.4 nextflow run metagenomics/metagenomics-tk \
	  -profile standard \
	  -params-file https://raw.githubusercontent.com/metagenomics/metagenomics-tk/refs/heads/master/default/tutorials/tutorial1/fullpipeline_annotation.yml \
	  -ansi-log false \
	  -entry wFullPipeline \
	  -resume \
	  --input.paired.path https://raw.githubusercontent.com/metagenomics/metagenomics-tk/refs/heads/master/test_data/tutorials/tutorial1/reads.tsv \
	  --logDir logs_annotation 

make check LOG_DIR=$(pwd)/logs_annotation
