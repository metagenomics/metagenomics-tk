# Run metabat
runMetaBat.sh !{contigs} !{bam}

# Create temporary directory
TEMP_DIR=$(mktemp -d -p .)

BIN_CONTIG_MAPPING=${TEMP_DIR}/bin_contig_mapping.tsv
for bin in $(basename !{contigs})*/bin* ; do
	BIN_NAME="!{sample}_$(basename ${bin})"

	# Get id of the bin (e.g get 2 of the bin SAMPLEID_bin.2.fa)
	ID=$(echo ${BIN_NAME} | rev | cut -d '.' -f 2 | rev)

	# Append bin id to every header
	seqkit replace  -p '(.*)' -r "\${1} MAG=${ID}" $bin > ${BIN_NAME}

	# Create old to new header table
	OLD_HEADERS=${TEMP_DIR}/old_headers.tsv
	NEW_HEADERS=${TEMP_DIR}/new_headers.tsv
	HEADERS_MAPPING=${TEMP_DIR}/headers_mapping.tsv
	grep ">" ${bin} | sed 's/>//g' > ${OLD_HEADERS}
	grep ">" ${BIN_NAME} | sed 's/>//g' > ${NEW_HEADERS}

	paste -d$'\t' ${OLD_HEADERS} ${NEW_HEADERS} >> ${HEADERS_MAPPING}

	# Create bin to contig mapping
	cat ${OLD_HEADERS} | sed "s/^/${BIN_NAME}\t/g" >> ${BIN_CONTIG_MAPPING}
done

# Add sample and bin information to every contig in depth file
BIN_DEPTH=${TEMP_DIR}/bin_depth.tsv
head -n 1 *.depth.txt  | sed "s/$/\tBIN_ID\tSAMPLE/g" > ${BIN_DEPTH}
join -t$'\t' -1 1 -2 2  <(tail -n +2  *.depth.txt | sort -k 1,1 ) <(sort -k 2,2 ${BIN_CONTIG_MAPPING}) \
			| sed "s/$/\t!{sample}/g" >> ${BIN_DEPTH}

# Compute bin statistics
SEQKIT_BIN_STATS=${TEMP_DIR}/seqkit_bin_stats.tsv
seqkit stat -Ta *.fa | sed '1s/^/SAMPLE\t/' | sed "1 ! s/^/!{sample}\t/" > ${SEQKIT_BIN_STATS}

# Add bin id to stats
SEQKIT_BIN_STATS_IDS=${TEMP_DIR}/seqkit_bin_stats_correct_ids.tsv
paste -d$'\t'  <(cut -f 2 ${SEQKIT_BIN_STATS} \
	        | sed 's/\.fa$//g' \
	        | tail -n +2 \
	        | sed '1s;^;BIN_ID\n;') ${SEQKIT_BIN_STATS} > ${SEQKIT_BIN_STATS_IDS}

# Compute average depth
BIN_AVG_DEPTH=${TEMP_DIR}/bin_avg_depth.tsv
for id in  $(cut -f 6 ${BIN_DEPTH} | tail -n +2 | sort | uniq); do  
        mean=$(awk -v bin=$id ' $6 == bin  { total += $3; count++ } END { print total/count }' ${BIN_DEPTH}); echo -e "$id\t$mean" ;
done | sed '1s/^/BIN_ID\tCOVERAGE\n/' > ${BIN_AVG_DEPTH}

# Add average depth to seqkit stats
join -t$'\t' --header -1 3 -2 1 ${SEQKIT_BIN_STATS_IDS} ${BIN_AVG_DEPTH} > !{sample}_bins_stats.tsv

# Replace old contig header with new one
BIN_DEPTH_CORRECT_IDS=!{sample}_contig_depth.tsv
head -n 1 *.depth.txt > ${BIN_DEPTH_CORRECT_IDS}
join -t$'\t' -1 1 -2 1 <(sort -k 1,1 ${HEADERS_MAPPING}) <(tail -n +2 *.depth.txt | sort -k 1,1) | cut -f 2- >> ${BIN_DEPTH_CORRECT_IDS}
