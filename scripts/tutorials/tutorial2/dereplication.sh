# This file will be referenced on the online wiki

#cd ~/mgcourse/

NXF_VER=25.04.2 nextflow run metagenomics/metagenomics-tk \
	  -profile standard \
	  -params-file https://raw.githubusercontent.com/metagenomics/metagenomics-tk/refs/heads/doc/exhibition-day/tutorial/default/tutorials/tutorial2/fullPipeline_dereplication.yml \
	  -ansi-log false \
	  -entry wAggregatePipeline \
	  -resume \
    --input output \
	  --logDir dereplication 

make check LOG_DIR=$(pwd)/dereplication
