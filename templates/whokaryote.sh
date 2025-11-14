gunzip -c !{gff} > unzipped.gff

# Whokaryote returns exit code 1 if no contig satisfies the required minimum contig length.
trap 'if [[ $? == 1 && ! -s output/contigs*.fasta ]]; then echo "Could not find contigs longer then required min length"; exit 0; fi' EXIT 

whokaryote.py --contigs !{fasta} --threads !{task.cpus} --gff unzipped.gff --outdir output !{params.steps.annotation.whokaryote.additionalParams} 

# Add sample and bin id column
addColumns() {
  local input_file="$1"
  local output_file="$2"

  if [ -s "$input_file" ]; then
    sed -e "1s/^/SAMPLE\tBIN_ID\t/" \
 	    -e "2,$ s/^/!{sample}\t/" \
	    -e "2,$ s/^/!{binID}\t/" $input_file > $output_file 
  fi 
}

# Add missing column header
sed -i '1i CONTIG' output/eukaryote_contig_headers.txt output/prokaryote_contig_headers.txt

# Output tsv instead of csv
sed -i 's/,/\t/g' output/featuretable.csv

EUKARYOTE_HEADER="!{binID}_eukaryote_contig_headers.tsv"

addColumns output/featuretable.csv !{binID}_featuretable.tsv
addColumns output/featuretable_predictions_T.tsv !{binID}_featuretable_predictions_T.tsv
addColumns output/eukaryote_contig_headers.txt ${EUKARYOTE_HEADER}
addColumns output/prokaryote_contig_headers.txt !{binID}_prokaryote_contig_headers.tsv

# In case Tiara is used
if [ -f output/tiara_pred.txt ]; then
	addColumns output/tiara_pred.txt !{binID}_tiara_pred.tsv
fi

IS_EUKARYOTE=true
if [ -s "${EUKARYOTE_HEADER}" ]; then
  IS_EUKARYOTE=true	
else
  IS_EUKARYOTE=false
fi
