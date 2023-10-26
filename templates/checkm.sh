
# Prepare checkm patch, output directory and output file name
mkdir out
FILE_ID=$(mktemp XXXXXXXX)
FILE=!{sample}_checkm_${FILE_ID}.tsv

# Check developer documentation
if [ -z "!{EXTRACTED_DB}" ]
then
  DATABASE=!{params.databases}/checkm
  LOCK_FILE=${DATABASE}/lock.txt

  echo '{"dataRoot": "!{params.databases}/checkm/out", "remoteManifestURL": "https://data.ace.uq.edu.au/public/CheckM_databases/", "manifestType": "CheckM", "remoteManifestName": ".dmanifest", "localManifestName": ".dmanifest"}' > /tmp/DATA_CONFIG

  if [ ! -z "!{S3_checkm_ACCESS}" ]
  then
    export AWS_ACCESS_KEY_ID=!{S3_checkm_ACCESS}
    export AWS_SECRET_ACCESS_KEY=!{S3_checkm_SECRET}
  fi

  # Download checkm database if necessary
  mkdir -p ${DATABASE}
  flock ${LOCK_FILE} concurrentDownload.sh --output=${DATABASE} \
    --link=!{DOWNLOAD_LINK} \
    --httpsCommand="wget -O checkm.tar.gz !{DOWNLOAD_LINK} && tar -xzvf checkm.tar.gz && rm checkm.tar.gz" \
    --s3FileCommand="s5cmd !{S5CMD_PARAMS} cp --concurrency !{task.cpus} !{DOWNLOAD_LINK} checkm.tar.gz && tar -xzvf checkm.tar.gz && rm checkm.tar.gz" \
    --s3DirectoryCommand="s5cmd !{S5CMD_PARAMS} cp --concurrency !{task.cpus} !{DOWNLOAD_LINK} . " \
    --s5cmdAdditionalParams="!{S5CMD_PARAMS}" \
    --localCommand="tar -xzvf !{DOWNLOAD_LINK}" \
    --expectedMD5SUM=!{MD5SUM}
else
  echo '{"dataRoot": "!{EXTRACTED_DB}", "remoteManifestURL": "https://data.ace.uq.edu.au/public/CheckM_databases/", "manifestType": "CheckM", "remoteManifestName": ".dmanifest", "localManifestName": ".dmanifest"}' > /tmp/DATA_CONFIG
fi

# run suggested checkm commands
checkm tree !{params.steps.magAttributes.checkm.additionalParams.tree} --pplacer_threads !{task.cpus} -t !{task.cpus} -x !{ending} . out > >(tee -a checkm_stdout.log) 2> >(tee -a checkm_stderr.log >&2)
checkm tree_qa out > >(tee -a checkm_stdout.log) 2> >(tee -a checkm_stderr.log >&2)
checkm lineage_set !{params.steps.magAttributes.checkm.additionalParams.lineage_set} out out/marker > >(tee -a checkm_stdout.log) 2> >(tee -a checkm_stderr.log >&2)
checkm analyze -x !{ending} -t !{task.cpus} out/marker . out > >(tee -a checkm_stdout.log) 2> >(tee -a checkm_stderr.log >&2)
checkm qa !{params.steps.magAttributes.checkm.additionalParams.qa} --tab_table -t !{task.cpus} -f checkm.txt out/marker out > >(tee -a checkm_stdout.log) 2> >(tee -a checkm_stderr.log >&2)

# If there is not enough disk space available, checkm exits with exit code 0
# The only solution is to check the stdout for errors.
if grep -q "Error" checkm_stdout.log checkm_stderr.log; then
	echo "Error found in CheckM log.";
	exit 1 ;
else
	echo "No error found in CheckM log.";
fi

# reformat output files according to magAttributes standard
echo -e "SAMPLE\tBIN_ID\tMarker lineage\t# genomes\t# markers\t# marker sets\t0\t1\t2\t3\t4\t5+\tCOMPLETENESS\tCONTAMINATION\tHETEROGENEITY" > $FILE
sed -i " 2,$ s/\t/!{ending}\t/"  checkm.txt
tail -n +2 checkm.txt | sed "s/^/!{sample}\t/g"  >> $FILE

