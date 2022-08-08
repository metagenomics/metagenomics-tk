for file in !{contigHeaderFiles}; do
	csvtk cut -f CONTIG --tabs ${file} | tail -n +2 >> filtered_tools_header.tsv
done

if [ -s filtered_tools_header.tsv ]; then
	sort filtered_tools_header.tsv <(seqkit fx2tab --name --only-id !{contigs}) \
	  | uniq -c | sed "s/^\ *//g" \
	  | grep "^!{NUMBER_OF_CONTIGS} " \
	  | cut -d ' ' -f 2- > filtered_selected_header.tsv


	if [ -s filtered_selected_header.tsv ]; then
		seqkit grep -f filtered_selected_header.tsv !{contigs} \
		 | seqkit seq --min-len !{MIN_LENGTH} \
		 | pigz -c > !{sample}_!{binID}_plasmids_filtered.fasta.gz 
	fi
fi
