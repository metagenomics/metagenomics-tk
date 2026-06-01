# --8<-- [start:func]
NXF_VER=25.10.4 nextflow run metagenomics/metagenomics-tk \
	  -profile standard \
	  -params-file https://raw.githubusercontent.com/metagenomics/metagenomics-tk/refs/heads/master/default/modules/binning/multiBinningONT.yml \
	  -ansi-log false \
	  -entry wMultiBinningLongRead \
	  -resume \
	  --steps.multiBinningONT.input.reads https://raw.githubusercontent.com/metagenomics/metagenomics-tk/refs/heads/master/test_data/multiBinningONT/samplesONT.tsv \
	  --steps.multiBinningONT.input.contigs https://raw.githubusercontent.com/metagenomics/metagenomics-tk/refs/heads/master/test_data/multiBinningONT/assemblyONT.tsv \
	  --steps.multiBinningONT.input.quality https://raw.githubusercontent.com/metagenomics/metagenomics-tk/refs/heads/master/test_data/multiBinningONT/quality.tsv \
	  --steps.multiBinningONT.input.groups https://raw.githubusercontent.com/metagenomics/metagenomics-tk/refs/heads/master/test_data/multiBinningONT/groups.tsv \
	  --logDir logs_multiBinningONT 
# --8<-- [end:func]

make check LOG_DIR=$(pwd)/logs_multiBinningONT
