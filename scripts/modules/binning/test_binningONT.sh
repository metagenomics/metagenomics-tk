# This file will be referenced on the online wiki

NXF_VER=25.10.4 nextflow run metagenomics/metagenomics-tk \
	  -profile standard \
	  -params-file https://raw.githubusercontent.com/metagenomics/metagenomics-tk/refs/heads/master/default/modules/binning/binningONT.yml \
	  -ansi-log false \
	  -entry wOntBinning \
	  -resume \
          --steps.binningONT.input.reads https://raw.githubusercontent.com/metagenomics/metagenomics-tk/refs/heads/master/test_data/binningONT/samplesONT.tsv \
          --steps.binningONT.input.contigs https://raw.githubusercontent.com/metagenomics/metagenomics-tk/refs/heads/master/test_data/binningONT/assemblyONT.tsv \
          --steps.binningONT.input.quality https://raw.githubusercontent.com/metagenomics/metagenomics-tk/refs/heads/master/test_data/binningONT/quality.tsv \
          --steps.binningONT.input.graph https://raw.githubusercontent.com/metagenomics/metagenomics-tk/refs/heads/master/test_data/binningONT/graph.tsv \
          --steps.binningONT.input.headerMapping https://raw.githubusercontent.com/metagenomics/metagenomics-tk/refs/heads/master/test_data/binningONT/headerMapping.tsv \
          --steps.binningONT.input.assemblyInfo https://raw.githubusercontent.com/metagenomics/metagenomics-tk/refs/heads/master/test_data/binningONT/headerMapping.tsv \
	  --logDir logs_binningONT 

make check LOG_DIR=$(pwd)/logs_binningONT
