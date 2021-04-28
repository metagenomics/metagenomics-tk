#!/bin/bash


# Check for Java
if ! command -v java &> /dev/null
then
	echo "Java is not installed"
	exit
fi

# Check for samtools
if ! command -v samtools &> /dev/null
then
	echo "Samtools is not installed"
	exit
fi


wget https://data.cami-challenge.org/camiClient.jar

java -jar camiClient.jar -d https://openstack.cebitec.uni-bielefeld.de:8080/swift/v1/CAMI_I_LOW data -p fq.gz


# prepare Data
test_binning_assembly_file=test_samples_interleaved.tsv
touch test_samples_interleaved.tsv
mkdir -p data
path=$(readlink -f .)
echo -e "SAMPLE\tREADS"  >> $test_binning_assembly_file
echo -e "test\t${path}/data/RL_S001__insert_270.fq.gz"  >> $test_binning_assembly_file

