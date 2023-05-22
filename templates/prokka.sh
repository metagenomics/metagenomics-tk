# Prepare Input Variables
BIN=!{fasta}
BIN_PREFIX=$(echo "${BIN%.*}")
BIN_ID="!{binID}"

# Run Prokka
if [[ !{fasta} == *.gz ]]; then
   zcat -f !{fasta} > input.fasta
   prokka !{params.steps.annotation.prokka.additionalParams} !{prodigalModeStr} --partialgenes --cpus !{task.cpus} --outdir out !{prokkaDomain} input.fasta
   rm input.fasta
else 
   prokka !{params.steps.annotation.prokka.additionalParams} !{prodigalModeStr} --partialgenes --cpus !{task.cpus} --outdir out !{prokkaDomain} !{fasta}
fi

# Rename files
for f in out/* ; do suffix=$(echo "${f##*.}"); mv $f ${BIN_PREFIX}.${suffix}; done

# create tsv out of gff
gffread ${BIN_PREFIX}.gff --table @id,@geneid,@chr,@start,@end,@strand,@numexons,@exons,@cds,@covlen,@cdslen \
     | sed '1 i\ID\tGENE_ID\tCHR\tSTART\tEND\tSTRAND\tNUMEXONS\tEXONS\tCDS\tCOVLEN\tCDSLEN' > gff.tsv

# join prokka tsv output with the tsv version of the gff file for additional contig information.
PROKKA_TMP_TSV=!{sample}_${BIN_ID}_prokka_tmp.tsv
PROKKA_COV_TMP_TSV=!{sample}_${BIN_ID}_prokka_cov_tmp.tsv
PROKKA_TSV=${BIN_ID}_prokka.tsv
csvtk join -d$'\t' -T -f "locus_tag;ID" ${BIN_PREFIX}.tsv gff.tsv > ${PROKKA_TMP_TSV}

# join prokka tsv with contig coverage tsv if not empty
if [ -s !{metabatCoverage} ] && [ -s !{defaultCoverage} ]; then 
	csvtk join -d$'\t' -T -f "CHR;CONTIG_NAME" ${PROKKA_TMP_TSV} !{metabatCoverage} \
		| csvtk cut -d$'\t' -T -f -SAMPLE  > ${PROKKA_COV_TMP_TSV}
	csvtk join -d$'\t' -T -f "CONTIG;CHR" !{defaultCoverage} ${PROKKA_COV_TMP_TSV} > ${PROKKA_TSV}
	else
	cp ${PROKKA_TMP_TSV} ${PROKKA_TSV}
fi

# Create SAMPLE, BIN_ID and CONTIG ids file
csvtk -t cut -f SAMPLE,CONTIG,locus_tag $PROKKA_TSV > prokka_tags_tmp.tsv

# Generate a sed script with substitution commands for each line in the prokka_tmp.tsv file
while read line; do
    # Extract the new tag from the current line
    newTag=$(echo -e "$line" | awk -v prefix="$BIN_PREFIX" -F '\t' '{gsub($1"_", "", $2); print prefix"_"$2"_"$3}')
    # Extract the old tag from the current line
    oldTag=$(echo -e "$line" | cut -f3)
    # Add a substitution command to the sed script
    echo "s/$oldTag/$newTag/g"
done < <(tail -n +2 prokka_tags_tmp.tsv) > sed_script.sed

# Find all files with the specified extensions and apply the sed script to each file
find . -type f \( -name "*.faa" -o -name "*.ffn" -o -name "*.gbk" -o -name "*.gff" -o -name "*.sqn" -o -name "*.tbl" \) -exec sed -i -f sed_script.sed {} +

# compress files individually to avoid errors if optional files are missing
for file in *; do
  if [[ $file == *.gff || $file == *.faa || $file == *.fna || $file == *.ffn || $file == *.fsa || $file == *.gbk || $file == *.sqn || $file == *tbl ]]; then
    pigz --best --processes !{task.cpus} "$file"
  fi
done
