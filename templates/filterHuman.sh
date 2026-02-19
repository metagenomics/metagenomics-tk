DB_FILE=""
if [ -z "${EXTRACTED_DB}" ]
then
	DATABASE=${params.polished.databases}/human_filter_db
	LOCK_FILE=\${DATABASE}/lock.txt

    	# Check if access and secret keys are necessary for s5cmd
    	if [ ! -z "${S3_filter_ACCESS}" ]
    	then
       		export AWS_ACCESS_KEY_ID=${S3_filter_ACCESS}
       		export AWS_SECRET_ACCESS_KEY=${S3_filter_SECRET}
    	fi

    	mkdir -p \${DATABASE}
    	flock \${LOCK_FILE} concurrentDownload.sh --output=\${DATABASE} \\
			--link=${DOWNLOAD_LINK} \\
			--httpsCommand="wgetStatic --no-check-certificate -qO- ${DOWNLOAD_LINK} | pigz -fdc > human_filter.db " \\
		        --s3FileCommand="s5cmd ${S5CMD_PARAMS} cat --concurrency ${task.cpus} ${DOWNLOAD_LINK} | pigz -fdc > human_filter.db  " \\
			--s3DirectoryCommand="s5cmd ${S5CMD_PARAMS} cp --concurrency ${task.cpus} ${DOWNLOAD_LINK} . && mv * human_filter.db " \\
			--s5cmdAdditionalParams="${S5CMD_PARAMS}" \\
	                --localCommand="gunzip -c ${DOWNLOAD_LINK} > ./human_filter.db " \\
		        --expectedMD5SUM=${MD5SUM}

     	DB_FILE="\${DATABASE}/out/human_filter.db"
     else
	DB_FILE="${EXTRACTED_DB}"
fi

INTERLEAVED_TMP=interleaved_tmp.fq
UNPAIRED_TMP=unpaired_tmp.fq

# srub.sh only accepts decompressed files
pigz -dc ${interleavedReads} > \${INTERLEAVED_TMP}
pigz -dc ${unpairedReads} > \${UNPAIRED_TMP}

# scrub.sh stores intermediate results in TMP. We can change the target directory by specifying TMPDIR variable.
export TMPDIR=\$(pwd)

INTERLEAVED_OUT=${sample}_interleaved.filtered.fq.gz
INTERLEAVED_REMOVED=${sample}_interleaved.removed.fq.gz

UNPAIRED_OUT=${sample}_unpaired.filtered.fq.gz
UNPAIRED_REMOVED=${sample}_unpaired.removed.fq.gz

scrub.sh ${ADDITIONAL_PARAMS} -x -r -s -u \${INTERLEAVED_REMOVED} -d \${DB_FILE} -i \${INTERLEAVED_TMP} -p ${task.cpus} -o - | pigz > \${INTERLEAVED_OUT}
scrub.sh ${ADDITIONAL_PARAMS} -x -r -u \${UNPAIRED_REMOVED} -d \${DB_FILE} -i \${UNPAIRED_TMP} -p ${task.cpus} -o - | pigz > \${UNPAIRED_OUT}

paste -d\$'\\t' <(echo -e "STATE\\nBEFORE") <(echo -e "SAMPLE\\n${sample}") <(seqkit stats \${INTERLEAVED_TMP} --all -T) > ${sample}_interleaved_summary_before.tsv
paste -d\$'\\t' <(echo -e "STATE\\nBEFORE") <(echo -e "SAMPLE\\n${sample}") <(seqkit stats \${UNPAIRED_TMP} --all -T) > ${sample}_unpaired_summary_before.tsv

paste -d\$'\\t' <(echo -e "STATE\\nAFTER") <(echo -e "SAMPLE\\n${sample}") <(seqkit stats \${INTERLEAVED_OUT} --all -T) > ${sample}_interleaved_summary_after.tsv
paste -d\$'\\t' <(echo -e "STATE\\nAFTER") <(echo -e "SAMPLE\\n${sample}") <(seqkit stats \${UNPAIRED_OUT} --all -T) > ${sample}_unpaired_summary_after.tsv
