# This file will be referenced on the online wiki

NXF_VER=24.04.4 nextflow run metagenomics/metagenomics-tk \
	  -profile standard \
	  -entry wFullPipeline \
          -params-file  https://raw.githubusercontent.com/metagenomics/metagenomics-tk/refs/heads/master/default/quickstart.yml \
	  --logDir logs \
	  --s3SignIn false \
	  --scratch false \
	  --output output \
	  --databases $(pwd)/databases


make check LOG_DIR=$(pwd)/logs
