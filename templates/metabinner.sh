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
/usr/local/bin/scripts/jgi_summarize_bam_contig_depths --minContigLength $MIN_LENGTH \
	--noIntraDepthVariance --outputDepth $DEPTH !{bam}
cut -f -1,4- $DEPTH > $COVERAGE

# Create a kmer profile 
python /usr/local/bin/scripts/gen_kmer.py $CONTIGS $((MIN_LENGTH-1)) ${KMER_SIZE}

# fix shebang in all perl files
sed -i '1 s/^.*$/#!\/usr\/bin\/env perl/' /usr/local/bin/auxiliary/*.pl

# run metabinner
bash /usr/local/bin/run_metabinner.sh -a $CONTIGS -o $RESULT -d $COVERAGE -k $(pwd)/*.csv -p /usr/local/bin -t !{task.cpus}

# Create bins
METABINNER_RESULT=${RESULT}/metabinner_res/metabinner_result.tsv 
OLD_BINS=${TEMP_DIR}/old_bins
BIN_CONTIG_MAPPING=!{sample}_bin_contig_mapping.tsv
echo -e "BIN_ID\tCONTIG\tBINNER" > ${BIN_CONTIG_MAPPING}
mkdir ${OLD_BINS}
for binId in $(cut -f 2 $METABINNER_RESULT | sort | uniq); do

	# Create bins with the old header and the new header
	BIN_NAME="!{sample}_bin.${binId}.fa"
	OLD_BIN_NAME=${OLD_BINS}/"!{sample}_bin.${binId}.fa"
	grep -e "$(printf '\t')${binId}$" $METABINNER_RESULT | cut -f 1 >> ${binId}.tsv 
	seqkit grep --threads !{task.cpus} -f ${binId}.tsv $CONTIGS \
		| seqkit replace -p "\s.+" > ${OLD_BIN_NAME}
	seqkit replace  -p '(.*)' -r "\${1} MAG=${binId}" ${OLD_BIN_NAME} > ${BIN_NAME}

	# Create bin to contig mapping and add the used binning tool to each line
	grep ">" ${OLD_BIN_NAME} | sed 's/>//g' \
		| sed "s/^/${BIN_NAME}\t/g;s/$/\tmetabinner/" >> ${BIN_CONTIG_MAPPING}
done


# return not binned fasta files
BINNED_IDS=binned.tsv
NOT_BINNED=!{sample}_notBinned.fa
grep -h ">" *.fa | tr -d ">" > ${BINNED_IDS}
if [ -s ${BINNED_IDS} ]; then
       # Get all not binned Ids
        seqkit grep -vf ${BINNED_IDS} !{contigs} \
         | seqkit replace  -p '(.*)' -r "\${1} MAG=NotBinned" > ${NOT_BINNED}
else
        seqkit replace  -p '(.*)' -r "\${1} MAG=NotBinned" !{contigs} > ${NOT_BINNED}
fi

# Add not binned contigs to mapping and add the used binning tool to each line
if [ -s ${NOT_BINNED} ]; then
		  	seqkit seq ${NOT_BINNED} -n -i | sed "s/^/NotBinned\t/g;s/$/\tmetabinner/" >> ${BIN_CONTIG_MAPPING}
fi

# Quickfix
# Explanation: the metabinner biocontainer does contain perl scripts that use a hardcoded path to
# perl interpreter /usr/bin/perl which does not exist in the container. Thats why the container is 
# run as root and the permissions of all output files must be updated.
chown -R 1000:1000 .
