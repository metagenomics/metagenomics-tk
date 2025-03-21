# This file will be referenced on the online wiki

NXF_VER=24.10.4 nextflow run metagenomics/metagenomics-tk \
	  -profile standard \
	  -params-file https://raw.githubusercontent.com/metagenomics/metagenomics-tk/refs/heads/master/default/tutorials/tutorial1/fullpipeline_binning.yml \
	  -ansi-log false \
	  -entry wFullPipeline \
	  -resume \
	  --logDir logs_binning \
	  --s3SignIn false \
	  --scratch false

make check LOG_DIR=$(pwd)/logs_binning
