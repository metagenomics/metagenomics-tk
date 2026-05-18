# This file will be referenced on the online wiki

# --8<-- [start:func]
NXF_VER=25.10.4 nextflow run metagenomics/metagenomics-tk \
	  -profile standard \
	  -params-file https://raw.githubusercontent.com/metagenomics/metagenomics-tk/refs/heads/master/default/modules/binning/binning.yml \
	  -ansi-log false \
	  -entry wShortReadBinning \
	  -resume \
	  --steps.binning.input.paired https://raw.githubusercontent.com/metagenomics/metagenomics-tk/refs/heads/master/test_data/binning/samples.tsv \
	  --steps.binning.input.single https://raw.githubusercontent.com/metagenomics/metagenomics-tk/refs/heads/master/test_data/binning/samplesUnpaired.tsv \
	  --steps.binning.input.contigs https://raw.githubusercontent.com/metagenomics/metagenomics-tk/refs/heads/master/test_data/binning/assembly.tsv \
	  --logDir logs_binning 
# --8<-- [end:func]

make check LOG_DIR=$(pwd)/logs_binning
