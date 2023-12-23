shopt -s nullglob

# Run metabat
jgi_summarize_bam_contig_depths --outputDepth coverage.txt *.bam
metabat2 -i !{contigs} -a coverage.txt -o bins/bin -t !{task.cpus} !{params.steps.binningHybrid.metabat.additionalParams}

# Create temporary directory
TEMP_DIR=$(mktemp -d -p .)

BIN_CONTIG_MAPPING=!{sample}_bin_contig_mapping.tsv
echo -e "BIN_ID\tCONTIG" > ${BIN_CONTIG_MAPPING}
for bin in $(find bins -name "bin*.fa"); do
	BIN_NAME="!{sample}_$(basename ${bin})"

	# Get id of the bin (e.g get 2 of the bin SAMPLEID_bin.2.fa)
	ID=$(echo ${BIN_NAME} | rev | cut -d '.' -f 2 | rev)

	# Append bin id to every header
	seqkit replace  -p '(.*)' -r "\${1} MAG=${ID}" $bin > ${BIN_NAME}

	# Create bin to contig mapping
	grep ">" ${bin} | sed 's/>//g' \
		| sed "s/^/${BIN_NAME}\t/g" >> ${BIN_CONTIG_MAPPING}
done

# return not binned fasta files
BINNED_IDS=binned.tsv
grep -h ">" bins/bin* | tr -d ">" > ${BINNED_IDS}
if [ -s ${BINNED_IDS} ]; then
        # Get all not binned Ids
        NOT_BINNED=!{sample}_notBinned.fa
	seqkit grep -vf ${BINNED_IDS} !{contigs} \
	 | seqkit replace  -p '(.*)' -r "\${1} MAG=NotBinned" > ${NOT_BINNED}
fi

# return not binned fasta files
BINNED_IDS=binned.tsv
NOT_BINNED=!{sample}_notBinned.fa
grep -h ">" bins/bin* | tr -d ">" > ${BINNED_IDS}
if [ -s ${BINNED_IDS} ]; then
	# Get all not binned Ids
	seqkit grep -vf ${BINNED_IDS} !{contigs} \
		| seqkit replace  -p '(.*)' -r "\${1} MAG=NotBinned" > ${NOT_BINNED}
else
	seqkit replace  -p '(.*)' -r "\${1} MAG=NotBinned" !{contigs} > ${NOT_BINNED}
fi

