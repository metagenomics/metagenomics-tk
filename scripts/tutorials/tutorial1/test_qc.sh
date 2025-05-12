# This file will be referenced on the online wiki

cd ~/mgcourse/

NXF_VER=25.04.0 nextflow run metagenomics/metagenomics-tk \
	  -profile standard \
	  -params-file https://raw.githubusercontent.com/metagenomics/metagenomics-tk/refs/heads/master/default/tutorials/tutorial1/fullpipeline_qc.yml \
	  -ansi-log false \
	  -entry wFullPipeline \
	  --input.paired.path https://raw.githubusercontent.com/metagenomics/metagenomics-tk/refs/heads/master/test_data/tutorials/tutorial1/reads.tsv \
	  --logDir logs_qc 

make check LOG_DIR=$(pwd)/logs_qc
