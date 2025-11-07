#!/bin/bash
# Create a script that creates for every bin a seperate file 

bin="$1"
FILE_NAME="${bin}_${2}"
csvtk -T -t filter2 -f "\$BIN_ID=='${bin}'" tmp/header_bin_id.tsv  \
 | csvtk -T -t cut -f CONTIG_NEW \
 | seqkit grep -f - tmp/all_updated_header.faa  >> ${FILE_NAME} ;
