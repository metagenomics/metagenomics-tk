# This file will be referenced on the online wiki

NXF_VER=24.04.0 nextflow run metagenomics/metagenomics-tk \
	  -profile standard \
          -params-file https://raw.githubusercontent.com/metagenomics/metagenomics-tk/refs/heads/master/default/tutorials/tutorial1/fullpipeline_qc.yml \
	  -ansi-log false \
	  -entry wFullPipeline \
	  --logDir logs_qc \
	  --s3SignIn false \
	  --scratch false \
	  --output output \
          --input.paired.path https://raw.githubusercontent.com/metagenomics/metagenomics-tk/refs/heads/master/test_data/tutorials/tutorial1/reads.tsv

make check LOG_DIR=$(pwd)/logs_qc
