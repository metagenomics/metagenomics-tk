shopt -s nullglob

# Run metabat
!{percentIdentity}  runMetaBat.sh !{metabatParams} !{contigs} !{bam}

# Create temporary directory
TEMP_DIR=$(basename $(mktemp))
mkdir ${TEMP_DIR}

BIN_CONTIG_MAPPING=!{sample}_bin_contig_mapping.tsv
echo -e "BIN_ID\tCONTIG\tBINNER" > ${BIN_CONTIG_MAPPING}
for bin in $(find $(basename !{contigs})* -name "bin*.fa"); do
	BIN_NAME="!{sample}_$(basename ${bin})"

	# Get id of the bin (e.g get 2 of the bin SAMPLEID_bin.2.fa)
	ID=$(echo ${BIN_NAME} | rev | cut -d '.' -f 2 | rev)

	# Append bin id to every header
	seqkit replace  -p '(.*)' -r "\${1} MAG=${ID}" $bin > ${BIN_NAME}

	# Create bin to contig mapping and add the used binner to each line
	grep ">" ${bin} | sed 's/>//g' \
		| sed "s/^/${BIN_NAME}\t/g;s/$/\tmetabat/" >> ${BIN_CONTIG_MAPPING}
done

# return not binned fasta files
BINNED_IDS=binned.tsv
grep -h ">" $(basename !{contigs})*/bin* | tr -d ">" > ${BINNED_IDS}
if [ -s ${BINNED_IDS} ]; then
        # Get all not binned Ids
        NOT_BINNED=!{sample}_notBinned.fa
	seqkit grep -vf ${BINNED_IDS} !{contigs} \
	 | seqkit replace  -p '(.*)' -r "\${1} MAG=NotBinned" > ${NOT_BINNED}
fi

# return not binned fasta files
BINNED_IDS=binned.tsv
NOT_BINNED=!{sample}_notBinned.fa
grep -h ">" $(basename !{contigs})*/bin* | tr -d ">" > ${BINNED_IDS}
if [ -s ${BINNED_IDS} ]; then
	# Get all not binned Ids
	seqkit grep -vf ${BINNED_IDS} !{contigs} \
		| seqkit replace  -p '(.*)' -r "\${1} MAG=NotBinned" > ${NOT_BINNED}
else
	seqkit replace  -p '(.*)' -r "\${1} MAG=NotBinned" !{contigs} > ${NOT_BINNED}
fi

# Add not binned contigs to mapping and add the used binning tool to each line
if [ -s ${NOT_BINNED} ]; then
		  	seqkit seq ${NOT_BINNED} -n -i | sed "s/^/NotBinned\t/g;s/$/\tmetabat/" >> ${BIN_CONTIG_MAPPING}
fi

# Fix for ownership issue https://github.com/nextflow-io/nextflow/issues/4565
chmod a+rw -R ${TEMP_DIR}
