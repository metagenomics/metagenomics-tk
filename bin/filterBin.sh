#!/bin/bash
# Create a script that creates for every bin a seperate file 

bin=$1;
FILE_NAME="${bin}_${2}"
cat tmp/concat.tsv | csvtk -T -t filter2 -f "\$BIN_ID=='${bin}'" >> ${FILE_NAME} ;

