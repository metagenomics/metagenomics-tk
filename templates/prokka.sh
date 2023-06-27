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

# join prokka tsv with contig coverage tsv if not empty, if empty the AnnotationFile module was called directly
# and the default coverage file is used, as contig information is not available
if [ -s !{metabatCoverage} ] && [ -s !{defaultCoverage} ]; then 
	csvtk join -d$'\t' -T -f "CHR;CONTIG_NAME" ${PROKKA_TMP_TSV} !{metabatCoverage} \
		| csvtk cut -d$'\t' -T -f -SAMPLE  > ${PROKKA_COV_TMP_TSV}
	csvtk join -d$'\t' -T -f "CONTIG;CHR" !{defaultCoverage} ${PROKKA_COV_TMP_TSV} > ${PROKKA_TSV}
	# Create SAMPLE, BIN_ID and CONTIG ids file
  csvtk -t cut -f SAMPLE,CONTIG,locus_tag $PROKKA_TSV > prokka_tags_tmp.tsv
	else
	cp ${PROKKA_TMP_TSV} ${PROKKA_TSV}
	csvtk -t cut -f locus_tag $PROKKA_TSV > prokka_tags_tmp.tsv
	# A dumy c tag is used as a placeholer, as contig information is not available
	sed -i "s/^/!{sample}\tc\t/g" prokka_tags_tmp.tsv
	sed -i '1s/.*/SAMPLE\tCONTIG\tlocus_tag/' prokka_tags_tmp.tsv
fi

# This script reads a prokka_tags_tmp.tsv file and generates a dictionary with old tags as keys and new tags as values.
# It finds all files with the specified extensions and applies the substitutions to each file separately.
renameProkkaTags.py $BIN_PREFIX !{task.cpus}

# compress files individually to avoid errors if optional files are missing
for file in *; do
  if [[ $file == *.gff || $file == *.faa || $file == *.fna || $file == *.ffn || $file == *.fsa || $file == *.gbk || $file == *.sqn || $file == *tbl ]]; then
    pigz --best --processes !{task.cpus} "$file"
  fi
done
