
# Prepare checkm patch, output directory and output file name
mkdir out
FILE_ID=$(mktemp XXXXXXXX)
FILE=!{sample}_checkm_${FILE_ID}.tsv

EXTRACTED_DB="!{params.steps?.magAttributes?.checkm?.database?.extractedDBPath}"

# Check developer documentation
if [[ $EXTRACTED_DB == "null" ]] 
then
  DATABASE=!{params.databases}/checkm
  LOCK_FILE=${DATABASE}/checksum.txt
  DOWNLOAD_LINK=!{params.steps.magAttributes.checkm.database.download.source}
  MD5SUM=!{params.steps.magAttributes.checkm.database.download.md5sum}
  S5CMD_PARAMS=!{params.steps?.magAttributes?.checkm?.download?.s5cmd?.params}

  echo '{"dataRoot": "!{params.databases}/checkm/out", "remoteManifestURL": "https://data.ace.uq.edu.au/public/CheckM_databases/", "manifestType": "CheckM", "remoteManifestName": ".dmanifest", "localManifestName": ".dmanifest"}' > /tmp/DATA_CONFIG

  # Download checkm database if necessary
  mkdir -p ${DATABASE}
  flock ${LOCK_FILE} concurrentDownload.sh --output=${DATABASE} \
    --link=$DOWNLOAD_LINK \
    --httpsCommand="wget -O checkm.tar.gz ${DOWNLOAD_LINK} && tar -xzvf checkm.tar.gz && rm checkm.tar.gz" \
    --s3FileCommand="s5cmd ${S5CMD_PARAMS} cp ${DOWNLOAD_LINK} checkm.tar.gz && tar -xzvf checkm.tar.gz && rm checkm.tar.gz" \
    --s3DirectoryCommand="s5cmd ${S5CMD_PARAMS} cp ${DOWNLOAD_LINK} . " \
    --s5cmdAdditionalParams="${S5CMD_PARAMS}" \
    --localCommand="tar -xzvf ${DOWNLOAD_LINK}" \
    --expectedMD5SUM=${MD5SUM}
else
  echo '{"dataRoot": ${EXTRACTED_DB}, "remoteManifestURL": "https://data.ace.uq.edu.au/public/CheckM_databases/", "manifestType": "CheckM", "remoteManifestName": ".dmanifest", "localManifestName": ".dmanifest"}' > /tmp/DATA_CONFIG
fi

# run suggested checkm commands
checkm tree !{params.steps.magAttributes.checkm.additionalParams.tree} --pplacer_threads !{task.cpus}  -t !{task.cpus} -x !{ending} . out
checkm tree_qa out
checkm lineage_set !{params.steps.magAttributes.checkm.additionalParams.lineage_set} out out/marker
checkm analyze -x !{ending} -t !{task.cpus} out/marker . out
checkm qa !{params.steps.magAttributes.checkm.additionalParams.qa} --tab_table -t !{task.cpus} -f checkm.txt out/marker out

# reformat output files according to magAttributes standard
echo -e "SAMPLE\tBIN_ID\tMarker lineage\t# genomes\t# markers\t# marker sets\t0\t1\t2\t3\t4\t5+\tCOMPLETENESS\tCONTAMINATION\tHETEROGENEITY" > $FILE
sed -i " 2,$ s/\t/.fa\t/"  checkm.txt
tail -n +2 checkm.txt | sed "s/^/!{sample}\t/g"  >> $FILE

