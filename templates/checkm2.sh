
# Prepare checkm patch, output directory and output file name
mkdir out
FILE_ID=$(mktemp XXXXXXXX)
OUTPUT=!{sample}_checkm2_${FILE_ID}.tsv

# Check developer documentation
if [ -z "!{EXTRACTED_DB}" ]
then
  DATABASE=!{params.polished.databases}/checkm2
  LOCK_FILE=${DATABASE}/lock.txt

  # Download checkm database if necessary
  mkdir -p ${DATABASE}
  flock ${LOCK_FILE} concurrentDownload.sh --output=${DATABASE} \
    --link=!{DOWNLOAD_LINK} \
    --httpsCommand="wget -O checkm2.tar.gz !{DOWNLOAD_LINK} && tar -xzvf checkm2.tar.gz && rm checkm2.tar.gz" \
    --s3FileCommand="s5cmd !{S5CMD_PARAMS} cp --concurrency !{task.cpus} !{DOWNLOAD_LINK} checkm2.tar.gz && tar -xzvf checkm2.tar.gz && rm checkm2.tar.gz" \
    --s3DirectoryCommand="s5cmd !{S5CMD_PARAMS} cp --concurrency !{task.cpus} !{DOWNLOAD_LINK} . " \
    --s5cmdAdditionalParams="!{S5CMD_PARAMS}" \
    --localCommand="tar -xzvf !{DOWNLOAD_LINK}" \
    --expectedMD5SUM=!{MD5SUM}
 
  export CHECKM2DB=$(find $DATABASE -name "*.dmnd")
else
  export CHECKM2DB=$(find !{EXTRACTED_DB} -name "*.dmnd")
fi

checkm2 predict !{params.steps.magAttributes.checkm2.additionalParams} --threads !{task.cpus} --input !{bins} -o out
mv out/quality_report.tsv ${OUTPUT}

# Prepare output

# Rename column names to upper letters
sed -i -e " 1,1 s/Name/BIN_ID/" \
    -e " 1,1 s/Completeness/COMPLETENESS/" \
    -e " 1,1 s/Contamination/CONTAMINATION/" ${OUTPUT}

# Add missing file ending
sed -i " 2,$ s/\t/!{ending}\t/" ${OUTPUT}

# Add sample id
sed -i -e "1s/^/SAMPLE\t/" -e "2,$ s/^/!{sample}\t/" ${OUTPUT}

# Add Heterogeneity column to make it compatible with older checkm result files
sed -i -e " 1,1 s/$/\tHETEROGENEITY/" \
	-e " 2,$ s/$/\t0/g" ${OUTPUT}
