
mkdir header
mkdir tmp
mkdir missing
for file in !{contigHeaderFiles}; do
        TMP_FILE=tmp/${file}
        FINAL_FILE=header/${file}
        MISSING_FILE=missing/${file}
        METHOD="$(echo ${file} | rev | cut -d '_' -f 1 | rev | cut -d '.' -f 1)"
	csvtk cut -f CONTIG --tabs ${file} | sed -e "2,$ s/$/\tTRUE/g"  -e "1 s/$/\t${METHOD}/g" > ${TMP_FILE}
	csvtk cut -f CONTIG --tabs ${file} | tail -n +2 >> filtered_tools_header.tsv
        if [ $(wc -l < ${TMP_FILE}) -gt 1 ]; then
              mv ${TMP_FILE} ${FINAL_FILE}
        else
              mv ${TMP_FILE} ${MISSING_FILE}
        fi
done

PLASMID_OUT_FASTA=!{binID}_filtered.fasta.gz 
PLASMID_OUT_TSV=!{binID}_filtered.tsv

if [ -s filtered_tools_header.tsv ]; then
	sort filtered_tools_header.tsv <(seqkit fx2tab --name --only-id !{contigs}) \
	  | uniq -c | sed "s/^\ *//g" \
	  | grep "^!{NUMBER_OF_CONTIGS} " \
	  | cut -d ' ' -f 2- > filtered_selected_header.tsv

	if [ -s filtered_selected_header.tsv ]; then

	        csvtk -t join -f 1 <(echo "CONTIG"; seqkit fx2tab --name --only-id !{contigs}) header/*  -k --na FALSE > !{binID}_detection_tools.tsv

		for file in $(ls -1 missing/) ; do  
                   MISS_METHOD=$(csvtk cut -f -CONTIG --tabs $file);  
                   sed -i -e "2,$ s/$/\tFALSE/g" -e "1 s/$/\t${MISS_METHOD}/g" !{binID}_detection_tools.tsv
                done

		seqkit grep -f filtered_selected_header.tsv !{contigs} \
		 | seqkit seq --min-len !{MIN_LENGTH} \
		 | pigz -c > ${PLASMID_OUT_FASTA}

                seqkit fx2tab -H --length --only-id --gc --name ${PLASMID_OUT_FASTA} > ${PLASMID_OUT_TSV}
	fi
fi
