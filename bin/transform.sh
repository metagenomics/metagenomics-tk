#!/bin/bash -ue

SEQS_INPUT=$1 
SEQS_OUTPUT=$2
HEADER_PREFIX=$3
CORES=$4

# The following function modifies the assembly fasta headers according to the pattern: SAMPLEID_SEQUENCECOUNTER_SEQUENCEHASH
	 
# Create temporary directory
TEMP_DIR=$(mktemp -d -p .)

# Extract assembly header
OLD_FASTA_HEADERS=${TEMP_DIR}/old_headers.tsv
seqkit fx2tab $SEQS_INPUT | cut -f 1 |  sed 's/^>//g' > ${OLD_FASTA_HEADERS}

# While loop for extracting sequence hashes based on the first 5 characters of md5sum
SEQUENCE_HASHES=${TEMP_DIR}/sequence_hashes.tsv
while read -r seq; do 
	printf %s "$seq" | md5sum | cut -f1 -d' '; 
done < <(seqkit fx2tab $SEQS_INPUT | cut -f 2) | cut -c -6 > ${SEQUENCE_HASHES}

NEW_FASTA_HEADERS=${TEMP_DIR}/new_headers.tsv
cat ${SEQUENCE_HASHES} \
	| awk '{print NR,$0}' \
	| sed "s/ /_/" \
	| sed "s/^/${HEADER_PREFIX}_/g" > ${NEW_FASTA_HEADERS}

# Concatenate SAMPLEID_SEQUENCECOUNTER with SEQUENCEHASH strings
FASTA_HEADERS_MAPPING=${TEMP_DIR}/header_mapping.tsv
csvtk concat --out-tabs -H <(csvtk transpose ${OLD_FASTA_HEADERS}) <(csvtk transpose ${NEW_FASTA_HEADERS}) \
	| csvtk --tabs transpose > ${FASTA_HEADERS_MAPPING}

# Replace old with new fasta header
seqkit replace -p '(.*)'  -r '{kv}' -k ${FASTA_HEADERS_MAPPING} ${SEQS_INPUT} \
		| pigz --best --processes $CORES > ${SEQS_OUTPUT} 