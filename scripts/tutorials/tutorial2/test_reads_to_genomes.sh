# This file will be referenced on the online wiki

cd ~/mgcourse/

NXF_VER=25.04.2 nextflow run metagenomics/metagenomics-tk \
	  -profile standard \
	  -params-file https://raw.githubusercontent.com/metagenomics/metagenomics-tk/refs/heads/master/default/tutorials/tutorial2/fullPipeline_reads_to_genomes.yml \
	  -ansi-log false \
	  -entry wFullPipeline \
	  -resume \
	  --input.paired.path https://raw.githubusercontent.com/metagenomics/metagenomics-tk/refs/heads/master/test_data/tutorials/tutorial2/samples.tsv \
	  --logDir logs_reads_to_genomes 

make check LOG_DIR=$(pwd)/logs_reads_to_genomes
