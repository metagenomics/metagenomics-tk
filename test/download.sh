#!/bin/bash
TARGET_PATH=$1


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

java -jar camiClient.jar -d https://openstack.cebitec.uni-bielefeld.de:8080/swift/v1/CAMI_I_LOW $TARGET_PATH -p fq.gz
