# run megahit
megahit -t !{task.cpus} --12 interleaved.fq.gz  !{includeUnpairedReads} !{params.steps.assembly.megahit.additionalParams}

ASSEMBLY_OUTPUT=!{sample}_contigs.fa.gz

# The following function modifies the assembly fasta headers according to the pattern: SAMPLEID_SEQUENCECOUNTER_SEQUENCEHASH
function setFastaHeader(){
	 
	# Create temporary directory
	TEMP_DIR=$(mktemp -d -p .)

	ASSEMBLY=megahit_out/final.contigs.fa 

	# Extract assembly header
	OLD_FASTA_HEADERS=${TEMP_DIR}/old_headers.tsv
	zgrep ">" $ASSEMBLY |  sed 's/^>//g' > ${OLD_FASTA_HEADERS}

	# While loop for extracting sequence hashes based on the first 5 characters of md5sum
	SEQUENCE_HASHES=${TEMP_DIR}/sequence_hashes.tsv
	while read -r seq; do 
		printf %s "$seq" | md5sum | cut -f1 -d' '; 
	done < <(seqkit fx2tab $ASSEMBLY | cut -f 2) | cut -c -6 > ${SEQUENCE_HASHES}

	NEW_FASTA_HEADERS=${TEMP_DIR}/new_headers.tsv
	cat ${SEQUENCE_HASHES} \
		| nl -w1 \
		| sed "s/\t/_/" \
		| sed "s/^/!{sample}_/g" > ${NEW_FASTA_HEADERS}

	# Concatenate SAMPLEID_SEQUENCECOUNTER with SEQUENCEHASH strings
	FASTA_HEADERS_MAPPING=${TEMP_DIR}/header_mapping.tsv
	paste -d$'\t' ${OLD_FASTA_HEADERS} ${NEW_FASTA_HEADERS} > ${FASTA_HEADERS_MAPPING}

	# Replace old with new fasta header
	seqkit replace -p '(.*)'  -r '{kv}' -k ${FASTA_HEADERS_MAPPING} ${ASSEMBLY} \
			| pigz --best --processes !{task.cpus} > ${ASSEMBLY_OUTPUT} 
}

setFastaHeader

# get basic contig stats 
paste -d$'\t' <(echo -e "SAMPLE\n!{sample}") <(seqkit stat -Ta ${ASSEMBLY_OUTPUT}) > !{sample}_contigs_stats.tsv

