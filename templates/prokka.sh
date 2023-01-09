# Prepare Input Variables
BIN=!{fasta}
BIN_PREFIX=$(echo "${BIN%.*}")
BIN_ID="$(basename !{fasta})"

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
PROKKA_TSV=!{sample}_${BIN_ID}_prokka.tsv
csvtk join -d$'\t' -T -f "locus_tag;ID" ${BIN_PREFIX}.tsv gff.tsv > ${PROKKA_TMP_TSV}

# join prokka tsv with contig coverage tsv
csvtk join -d$'\t' -T -f "CHR;CONTIG_NAME" ${PROKKA_TMP_TSV} !{metabatCoverage} \
	| csvtk cut -d$'\t' -T -f -SAMPLE  > ${PROKKA_COV_TMP_TSV}
csvtk join -d$'\t' -T -f "CONTIG;CHR" !{defaultCoverage} ${PROKKA_COV_TMP_TSV} > ${PROKKA_TSV}

pigz --best --processes !{task.cpus} *gff *.faa *.fna *.ffn *.fsa *.gbk *.sqn *tbl
