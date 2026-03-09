mkdir output

SemiBin2 single_easy_bin --compression none --processes ${task.cpus} -i ${contigs} -b ${bam} -o output ${semibinParams}

BIN_CONTIG_MAPPING=${sample}_bin_contig_mapping.tsv
echo -e "BIN_ID\tCONTIG\tBINNER" > \${BIN_CONTIG_MAPPING}
for bin in \$(find output/output_bins/ -name "*.fa"); do 

	# Get id of the bin (e.g get 2 of the bin SemiBin_2.fa)
	ID=\$(echo \$(basename \$bin) | rev | cut -d '_' -f 1 | rev | cut -d '.' -f 1  )
	BIN_NAME="${sample}_bin.\${ID}.fa"
  
  # Append bin id to every header
	seqkit replace  -p '(.*)' -r "\\${1} MAG=\${ID}" \$bin > \${BIN_NAME}

	# Create bin to contig mapping and add the used binner to each line
	grep ">" \${bin} | sed 's/>//g' \\
		| sed "s/^/\${BIN_NAME}\t/g;s/\$/\tsemibin2/" >> \${BIN_CONTIG_MAPPING}
done

# return not binned fasta files
BINNED_IDS=binned.tsv
NOT_BINNED=${sample}_notBinned.fa
grep -h ">" \$(find output/output_bins/ -name "*.fa") | tr -d ">" > \${BINNED_IDS}
if [ -s \${BINNED_IDS} ]; then
	# Get all not binned Ids
	seqkit grep -vf \${BINNED_IDS} ${contigs} \\
		| seqkit replace  -p '(.*)' -r "\\${1} MAG=NotBinned" > \${NOT_BINNED}
else
	seqkit replace  -p '(.*)' -r "\\${1} MAG=NotBinned" ${contigs} > \${NOT_BINNED}
fi

# Fix for ownership issue https://github.com/nextflow-io/nextflow/issues/4565
chmod a+rw -R output
