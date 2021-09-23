DEPTH="$(basename !{contigs}).depth.txt"
COVERAGE="$(pwd)/coverage_profile.tsv"
metabinner_path=$(dirname $(which run_metabinner.sh))
CONTIGS="$(pwd)/contigs.fasta"
RESULT="$(pwd)/result"

MIN_LENGTH=!{params.steps.binning.metabinner.minContigLength}
KMER_SIZE=!{params.steps.binning.metabinner.kmerSize}

# Create temporary directory
TEMP_DIR=$(mktemp -d -p .)

# Metabinner works only with unzipped fasta files
gunzip -dc !{contigs} | seqkit seq --quiet -m ${MIN_LENGTH} > $CONTIGS

# Create coverage depth file based on the alignment
/usr/local/bin/scripts/jgi_summarize_bam_contig_depths --minContigLength $MIN_LENGTH --noIntraDepthVariance --outputDepth $DEPTH !{bam}
cut -f -1,4- $DEPTH > $COVERAGE

# Create a kmer profile 
python /usr/local/bin/scripts/gen_kmer.py $CONTIGS $((MIN_LENGTH-1)) ${KMER_SIZE}

# fix shebang in all perl files
sed -i '1 s/^.*$/#!\/usr\/bin\/env perl/' /usr/local/bin/auxiliary/*.pl

# run metabinner
bash /usr/local/bin/run_metabinner.sh -a $CONTIGS -o $RESULT -d $COVERAGE -k $(pwd)/*.csv -p /usr/local/bin -t !{task.cpus}

# Run jgi_summarize_bam_contig_depths again in order to make the coverage information comparable to other tools
rm $DEPTH
/usr/local/bin/scripts/jgi_summarize_bam_contig_depths --minContigLength $MIN_LENGTH --outputDepth $DEPTH !{bam}
cut -f -1,4- $DEPTH > $COVERAGE

# Create bins
METABINNER_RESULT=${RESULT}/metabinner_res/metabinner_result.tsv 
OLD_BINS=${TEMP_DIR}/old_bins
BIN_CONTIG_MAPPING=${TEMP_DIR}/bin_contig_mapping.tsv
HEADERS_MAPPING=${TEMP_DIR}/headers_mapping.tsv
mkdir ${OLD_BINS}
for binId in $(cut -f 2 $METABINNER_RESULT | sort | uniq); do

	# Create bins with the old header and the new header
	BIN_NAME="!{sample}_bin.${binId}.fa"
	OLD_BIN_NAME=${OLD_BINS}/"!{sample}_bin.${binId}.fa"
	grep -e "$(printf '\t')${binId}$" $METABINNER_RESULT | cut -f 1 >> ${binId}.tsv 
	seqkit grep -f ${binId}.tsv $CONTIGS \
		| seqkit replace -p "\s.+" > ${OLD_BIN_NAME}
	seqkit replace  -p '(.*)' -r "\${1}_${binId}" ${OLD_BIN_NAME} > ${BIN_NAME}

	# Create old to new header table
	OLD_HEADERS=${TEMP_DIR}/old_headers.tsv
	NEW_HEADERS=${TEMP_DIR}/new_headers.tsv
	grep ">" ${OLD_BIN_NAME} | sed 's/>//g' > ${OLD_HEADERS}
	grep ">" ${BIN_NAME} | sed 's/>//g' > ${NEW_HEADERS}

	paste -d$'\t' ${OLD_HEADERS} ${NEW_HEADERS} >> ${HEADERS_MAPPING}

	# Create bin to contig mapping
	cat ${OLD_HEADERS} | sed "s/^/${BIN_NAME}\t/g" >> ${BIN_CONTIG_MAPPING}
done

# Add sample and bin information to every contig in depth file
#BIN_DEPTH=${TEMP_DIR}/bin_depth.tsv

BIN_DEPTH=!{sample}_bins_depth.tsv
head -n 1 *.depth.txt  | sed "s/$/\tBIN_ID\tSAMPLE/g" > ${BIN_DEPTH}
csvtk join -t -H -f "1;2" <(tail -n +2 *.depth.txt) <(tail -n +2 ${BIN_CONTIG_MAPPING}) \
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
	        mean=$(awk -v bin=$id ' $6 == bin  { total += $3; count++ } END { print total/count }' ${BIN_DEPTH})
		echo -e "$id\t$mean" 
done | sed '1s/^/BIN_ID\tCOVERAGE\n/' > ${BIN_AVG_DEPTH}

# Add average depth to seqkit stats
csvtk  join -t -f "3;1" ${SEQKIT_BIN_STATS_IDS} ${BIN_AVG_DEPTH} > !{sample}_bins_stats.tsv

# Quickfix
# Explanation: the metabinner biocontainer does contain perl scripts that use a hardcoded path to
# perl interpreter /usr/bin/perl which does not exist in the container. Thats why the container is 
# run as root and the permissions of all output files must be updated.
chown -R 1000:1000 .
