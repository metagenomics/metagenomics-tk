# Create temporary directory
TEMP_DIR=$(mktemp -d -p .)

DEPTH=!{sample}.depth.tsv
jgi_summarize_bam_contig_depths !{bam} --outputDepth $DEPTH

# Add sample and bin information to every contig in depth file
BIN_CONTIG_MAPPING=!{binContigMapping}
BIN_DEPTH=!{sample}_contigs_depth.tsv
head -n 1 *.depth.tsv  | sed "s/$/\tBIN_ID\tSAMPLE/g" > ${BIN_DEPTH}
csvtk join -t -H -f "1;2" <(tail -n +2 *.depth.tsv) <(tail -n +2 ${BIN_CONTIG_MAPPING}) \
	        | sed "s/$/\t!{sample}/g" >> ${BIN_DEPTH}

# Compute bin statistics
SEQKIT_BIN_STATS=${TEMP_DIR}/seqkit_bin_stats.tsv
seqkit stat --threads !{task.cpus} -Ta !{bins} | sed '1s/^file/BIN_ID/' | sed '1s/^/SAMPLE\t/' \
	| sed "1 ! s/^/!{sample}\t/" > ${SEQKIT_BIN_STATS}

# Compute average depth
BIN_AVG_DEPTH=${TEMP_DIR}/bin_avg_depth.tsv
for id in $(cut -f 6 ${BIN_DEPTH} | tail -n +2 | sort | uniq); do
       mean=$(awk -v bin=$id ' $6 == bin  { total += $3; count++ } END { print total/count }' ${BIN_DEPTH})
       echo -e "$id\t$mean" 
done | sed '1s/^/BIN_ID\tCOVERAGE\n/' > ${BIN_AVG_DEPTH}

# Add average depth to seqkit stats
csvtk  join -t -f "BIN_ID" ${SEQKIT_BIN_STATS} ${BIN_AVG_DEPTH} > !{sample}_bins_stats.tsv
