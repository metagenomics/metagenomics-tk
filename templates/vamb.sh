VAMB_PARAMS="${params.steps.binning.vamb.additionalParams}"

# Create temporary directory
TEMP_DIR=\$(mktemp -d -p .)

vamb bin default --outdir vambout --fasta ${contigs} --bamdir . \${VAMB_PARAMS}
 
BIN_CONTIG_MAPPING=${sample}_bin_contig_mapping.tsv
echo -e "BIN_ID\tCONTIG\tBINNER" > \${BIN_CONTIG_MAPPING}
for bin in \$(find vambout/bins/ -name "*.fna"); do
 	BIN_NAME_OLD="\$(basename \${bin})"
 	ID=\$(echo \${BIN_NAME_OLD} | rev | cut -d '.' -f 2 | rev)
 	BIN_NAME="${sample}_bin.\${ID}.fa"

	# Append bin id to every header
 	seqkit replace  -p '(.*)' -r "\\\${1} MAG=\${ID}" \$bin > \${BIN_NAME}
 
 	# Create bin to contig mapping and add the used binner to each line
 	grep ">" \${bin} | sed 's/>//g' \\
 		| sed "s|^|\${BIN_NAME}\\t|g;s|\$|\\tVAMB|" >> \${BIN_CONTIG_MAPPING}
done
 
# return not binned fasta files
BINNED_IDS=binned.tsv
NOT_BINNED=${sample}_notBinned.fa
grep -h ">" \$(find vambout/bins/ -name "*.fna") | tr -d ">" > \${BINNED_IDS}
if [ -s \${BINNED_IDS} ]; then
	# Get all not binned Ids
	seqkit grep -vf \${BINNED_IDS} ${contigs} \\
		| seqkit replace  -p '(.*)' -r "\\\${1} MAG=NotBinned" > \${NOT_BINNED}
else
	seqkit replace  -p '(.*)' -r "\\\${1} MAG=NotBinned" ${contigs} > \${NOT_BINNED}
fi

# Fix for ownership issue https://github.com/nextflow-io/nextflow/issues/4565
chmod a+rw -R \${TEMP_DIR}
