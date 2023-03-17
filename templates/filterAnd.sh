
mkdir header
for file in !{contigHeaderFiles}; do
        METHOD="$(echo ${file} | rev | cut -d '_' -f 1 | rev | cut -d '.' -f 1)"
        csvtk cut -f CONTIG --tabs ${file} | sed -e "2,$ s/$/\tTRUE/g"  -e "1 s/$/\t${METHOD}/g" > header/${file}
	csvtk cut -f CONTIG --tabs ${file} | tail -n +2 >> filtered_tools_header.tsv
done

PLASMID_OUT_FASTA=!{binID}_filtered.fasta.gz 
PLASMID_OUT_TSV=!{binID}_filtered.tsv

if [ -s filtered_tools_header.tsv ]; then
	sort filtered_tools_header.tsv <(seqkit fx2tab --name --only-id !{contigs}) \
	  | uniq -c | sed "s/^\ *//g" \
	  | grep "^!{NUMBER_OF_CONTIGS} " \
	  | cut -d ' ' -f 2- > filtered_selected_header.tsv

	if [ -s filtered_selected_header.tsv ]; then

	        csvtk -t join -f 1 <(seqkit fx2tab --name --only-id !{contigs}) header/*  -k --na FALSE > !{binID}_detection_tools.tsv

		seqkit grep -f filtered_selected_header.tsv !{contigs} \
		 | seqkit seq --min-len !{MIN_LENGTH} \
		 | pigz -c > ${PLASMID_OUT_FASTA}

                seqkit fx2tab -H --length --only-id --gc --name ${PLASMID_OUT_FASTA} > ${PLASMID_OUT_TSV}
	fi
fi
