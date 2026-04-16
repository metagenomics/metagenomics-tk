# This file will be referenced on the online wiki

NXF_VER=25.10.4 nextflow run metagenomics/metagenomics-tk \
	  -profile standard \
	  -entry wFullPipeline \
	  -params-file  https://raw.githubusercontent.com/metagenomics/metagenomics-tk/refs/heads/master/default/quickstart.yml \
	  --logDir logs \
	  --s3SignIn false \
	  --scratch false \
	  --output output \
	  --databases $(pwd)/databases \
	  --input.paired.r1 https://openstack.cebitec.uni-bielefeld.de:8080/swift/v1/meta_test/small/read1_1.fq.gz \
	  --input.paired.r2 https://openstack.cebitec.uni-bielefeld.de:8080/swift/v1/meta_test/small/read2_1.fq.gz \
	  --input.paired.names test1

make check LOG_DIR=$(pwd)/logs
