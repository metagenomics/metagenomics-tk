#!/bin/bash
# Create a script that creates for every bin a seperate file 

bin=$1;
FILE_NAME="${bin}_${2}"
FILE_TO_CONCAT=$3
cat ${FILE_TO_CONCAT} | csvtk -T -t filter2 -f "\$BIN_ID=='${bin}'" >> ${FILE_NAME} ;
