# Check developer documentation
# Download titles
if [ -z "!{TITLES_EXTRACTED_DB}" ]
then
  TITLES_DATABASE=!{params.databases}/emgb_titles
  TITLES_LOCK_FILE=${TITLES_DATABASE}/lock.txt

  if [ ! -z "!{S3_EMGB_TITLES_ACCESS}" ]
  then
    export AWS_ACCESS_KEY_ID=!{S3_EMGB_TITLES_ACCESS}
    export AWS_SECRET_ACCESS_KEY=!{S3_EMGB_TITLES_SECRET}
  fi

  echo "${TITLES_DATABASE}"
  # Download checkm database if necessary
  mkdir -p ${TITLES_DATABASE}
  flock ${TITLES_LOCK_FILE} concurrentDownload.sh --output=${TITLES_DATABASE} \
    --link="!{TITLES_DOWNLOAD_LINK}" \
    --httpsCommand="wgetStatic --no-check-certificate -qO- !{TITLES_DOWNLOAD_LINK} | gunzip > titles.tsv" \
    --s3FileCommand="s5cmd !{TITLES_S5CMD_PARAMS} cat --concurrency !{task.cpus} !{TITLES_DOWNLOAD_LINK} | gunzip > titles.tsv " \
    --s3DirectoryCommand="s5cmd !{TITLES_S5CMD_PARAMS} cp --concurrency !{task.cpus} !{TITLES_DOWNLOAD_LINK} . && mv * titles.tsv " \
    --s5cmdAdditionalParams="!{TITLES_S5CMD_PARAMS}" \
    --localCommand="gunzip -c !{TITLES_DOWNLOAD_LINK} > ./titles.tsv" \
    --expectedMD5SUM=!{TITLES_MD5SUM}
  TITLES_FILE="${TITLES_DATABASE}/out/titles.tsv"
else
  TITLES_FILE="!{TITLES_EXTRACTED_DB}"
fi



# Check developer documentation
# Download titles
if [ -z "!{KEGG_EXTRACTED_DB}" ]
then
  KEGG_DATABASE=!{params.databases}/emgb_kegg
  KEGG_LOCK_FILE=${KEGG_DATABASE}/lock.txt

  if [ ! -z "!{S3_EMGB_KEGG_ACCESS}" ]
  then
    export AWS_ACCESS_KEY_ID=!{S3_EMGB_KEGG_ACCESS}
    export AWS_SECRET_ACCESS_KEY=!{S3_EMGB_KEGG_SECRET}
  fi

  # Download checkm database if necessary
  mkdir -p ${KEGG_DATABASE}
  flock ${KEGG_LOCK_FILE} concurrentDownload.sh --output=${KEGG_DATABASE} \
    --link=!{KEGG_DOWNLOAD_LINK} \
    --httpsCommand="wgetStatic --no-check-certificate -qO- !{KEGG_DOWNLOAD_LINK} | zstd -T!{task.cpus} -d -c | tar --strip-components=1 -xv " \
    --s3FileCommand="s5cmd !{KEGG_S5CMD_PARAMS} cat --concurrency !{task.cpus} !{KEGG_DOWNLOAD_LINK} | zstd -T!{task.cpus} -d -c | tar --strip-components=1 -xv " \
    --s3DirectoryCommand="s5cmd !{KEGG_S5CMD_PARAMS} cp --concurrency !{task.cpus} !{KEGG_DOWNLOAD_LINK} . " \
    --s5cmdAdditionalParams="!{KEGG_S5CMD_PARAMS}" \
    --localCommand="zstd -T!{task.cpus} -c -d !{KEGG_DOWNLOAD_LINK} | tar --strip-components=1 -xv " \
    --expectedMD5SUM=!{KEGG_MD5SUM}
  KEGG_DB=$(readlink -f ${KEGG_DATABASE}/out/)
else
  KEGG_DB="!{KEGG_EXTRACTED_DB}"
fi

BINS_DIR=bins

nr=$(find  blastResult/ -name "*.tsv" -exec readlink -f {} \;  | sed 's/^/ --blast-tab /g')

tax=$(find taxonomy/ -name "*.tsv" -exec readlink -f {} \; | sed 's/^/ -mmseqs-lineage /g')


ffn=$(find ffn/ -name "*.gz" -exec readlink -f {} \; | sed 's/^/ -ffn /g')

faa=$(find faa/ -name "*.gz" -exec readlink -f {} \; | sed 's/^/ -faa /g')

gff=$(find gff/ -name "*.gz" -exec readlink -f {} \; | sed 's/^/ -gff /g')

kegg=$(find blastKeggResult/ -name "*.tsv" -exec readlink -f {} \; | sed 's/^/ -kegg-blast-tab /g')

titles=" -title-tsv $TITLES_FILE "

name=" -dataset-name !{sample} "

db=" -db ${KEGG_DB} "

annotatedgenes2json -bins-dir bins \
			${nr} \
			${tax} \
			${ffn} \
			${faa} \
			${gff} \
			${titles} \
			${name} \
			${db} \
			${kegg} \
			-json-gz !{sample}.genes.json.gz
