# This file will be referenced on the online wiki

NXF_VER=25.10.4 nextflow run metagenomics/metagenomics-tk \
	  -profile standard \
	  -params-file https://raw.githubusercontent.com/metagenomics/metagenomics-tk/refs/heads/master/default/modules/binning/multiBinning.yml \
	  -ansi-log false \
	  -entry wMultiBinningShortRead \
	  -resume \
          --steps.multiBinning.input.paired https://raw.githubusercontent.com/metagenomics/metagenomics-tk/refs/heads/master/test_data/multiBinning/samples.tsv \
          --steps.multiBinning.input.single https://raw.githubusercontent.com/metagenomics/metagenomics-tk/refs/heads/master/test_data/multiBinning/samplesUnpaired.tsv \
          --steps.multiBinning.input.contigs https://raw.githubusercontent.com/metagenomics/metagenomics-tk/refs/heads/master/test_data/multiBinning/assembly.tsv \
          --steps.multiBinning.input.groups https://raw.githubusercontent.com/metagenomics/metagenomics-tk/refs/heads/master/test_data/multiBinning/groups.tsv \
	  --logDir logs_multiBinning 

make check LOG_DIR=$(pwd)/logs_multiBinning
