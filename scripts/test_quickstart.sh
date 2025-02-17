# This file will be referenced on the online wiki

NXF_VER=23.10.0 nextflow run metagenomics/metagenomics-tk \
	  -profile standard \
	  -entry wFullPipeline \
	  -params-file default/quickstart.yml \
	  --logDir logs \
	  --s3SignIn false \
	  --scratch false \
	  --output output \
	  --databases $(pwd)/databases


make check LOG_DIR=$(pwd)/logs
