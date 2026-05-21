# --8<-- [start:func]
NXF_VER=25.10.4 nextflow run metagenomics/metagenomics-tk \
	  -profile standard \
	  -params-file https://raw.githubusercontent.com/metagenomics/metagenomics-tk/refs/heads/master/default/modules/binning/binningONT.yml \
	  -ansi-log false \
	  -entry wOntBinning \
	  -resume \
	  --steps.binningONT.input.reads https://raw.githubusercontent.com/metagenomics/metagenomics-tk/refs/heads/master/test_data/binningONT/samplesONT.tsv \
	  --steps.binningONT.input.contigs https://raw.githubusercontent.com/metagenomics/metagenomics-tk/refs/heads/master/test_data/binningONT/assemblyONT.tsv \
	  --steps.binningONT.input.quality https://raw.githubusercontent.com/metagenomics/metagenomics-tk/refs/heads/master/test_data/binningONT/quality.tsv \
	  --logDir logs_binningONT 
# --8<-- [end:func]

make check LOG_DIR=$(pwd)/logs_binningONT
