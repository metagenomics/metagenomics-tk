# if no local database is referenced, start download part
if [ -z "!{EXTRACTED_DB}" ]; then
  # polished params end with a / so no additional one is needed at the end
  DATABASE=!{params.polished.databases}metabuli
  mkdir -p ${DATABASE}/!{dbType}
  LOCK_FILE=${DATABASE}/!{dbType}/lock.txt

  # Check if access and secret keys are necessary for s5cmd
  if [ ! -z "!{S3_TAX_DB_ACCESS}" ]; then
    export AWS_ACCESS_KEY_ID=!{S3_TAX_DB_ACCESS}
    export AWS_SECRET_ACCESS_KEY=!{S3_TAX_DB_SECRET}
  fi

  flock ${LOCK_FILE} concurrentDownload.sh --output=${DATABASE}/!{dbType} \
    --link="!{DOWNLOAD_LINK}" \
    --httpsCommand="wgetStatic --no-check-certificate -qO- !{DOWNLOAD_LINK} | tar -xzv " \
    --s3FileCommand="s5cmd !{S5CMD_PARAMS} cat --concurrency !{task.cpus} !{DOWNLOAD_LINK} | tar -xzv " \
    --s3DirectoryCommand="s5cmd !{S5CMD_PARAMS} cp --concurrency !{task.cpus} !{DOWNLOAD_LINK} . " \
    --s5cmdAdditionalParams="!{S5CMD_PARAMS}" \
    --localCommand="tar -xzvf !{DOWNLOAD_LINK}" \
    --expectedMD5SUM=!{MD5SUM}

  # The archive contains one directory, the Metabuli database. Its db.parameters file marks the database root.
  METABULI_DATABASE_DIR=$(dirname "$(find ${DATABASE}/!{dbType}/out -name db.parameters -type f | head -n 1)")
else
  METABULI_DATABASE_DIR="!{EXTRACTED_DB}"
fi

# Fix for ownership issue https://github.com/nextflow-io/nextflow/issues/4565
mkdir -p out tmp
chmod a+rw -R out tmp

metabuli classify --seq-mode 3 --lineage 1 \
  "!{contigs}" "${METABULI_DATABASE_DIR}" out "!{sample}_!{binType}_!{binId}.!{dbType}" !{parameters} \
  --threads !{task.cpus} \
  --max-ram $(echo !{task.memory} | cut -d ' ' -f 1)

# Add SAMPLE and BIN_ID as first columns
sed -i -e '1 s/^/SAMPLE\tBIN_ID\t/' -e "2,\$ s/^/!{sample}\t!{binId}\t/" out/*_classifications.tsv
