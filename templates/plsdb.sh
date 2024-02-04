
EXTRACTED_DB="!{params.steps?.plasmid?.PLSDB?.database?.extractedDBPath}"

# Check developer documentation
PLSDB=""
if [[ $EXTRACTED_DB == "null" ]]
then 
  DATABASE=!{params.databases}/plsdb
  LOCK_FILE=${DATABASE}/checksum.txt
  DOWNLOAD_LINK=!{params?.steps?.plasmid?.PLSDB?.database?.download?.source}
  MD5SUM=!{params?.steps?.plasmid?.PLSDB?.database?.download?.md5sum}

  # Check if access and secret keys are necessary for s5cmd
  if [ ! -z "!{S3_PLSDB_ACCESS}" ]
  then
    export AWS_ACCESS_KEY_ID=!{S3_PLSDB_ACCESS}
    export AWS_SECRET_ACCESS_KEY=!{S3_PLSDB_SECRET}
  fi

  # Download plsdb database if necessary
  mkdir -p ${DATABASE}
  flock ${LOCK_FILE} concurrentDownload.sh --output=${DATABASE} \
	--link=$DOWNLOAD_LINK \
	--httpsCommand="wget -qO- $DOWNLOAD_LINK | tar xjv " \
	--s3FileCommand="s5cmd !{S5CMD_PARAMS} cat --concurrency !{task.cpus} ${DOWNLOAD_LINK} |  tar xjv " \
	--s3DirectoryCommand="s5cmd !{S5CMD_PARAMS} cp --concurrency !{task.cpus} ${DOWNLOAD_LINK} . " \
	--s5cmdAdditionalParams="!{S5CMD_PARAMS}" \
	--localCommand="tar xjvf ${DOWNLOAD_LINK} " \
	--expectedMD5SUM=${MD5SUM}

   PLSDB=${DATABASE}/out
else
   PLSDB=${EXTRACTED_DB}
fi

# Sketch plasmids
mash sketch -i !{params.steps.plasmid.PLSDB.additionalParams.mashSketch} -o query !{plasmids}

# Compute distance between the query and plsdb plasmids
mash dist ${PLSDB}/plsdb.msh query.msh -p !{task.cpus} !{params.steps.plasmid.PLSDB.additionalParams.mashDist} > !{binID}.tsv

# filter matches by user defined threshold
MATCHES_RAW=matches.tsv
sort -rgk 5,5 !{binID}.tsv | sed 's/\/.*$//g' \
     | awk -v threshold=!{params.steps.plasmid.PLSDB.sharedKmerThreshold} '($5+0 > threshold) {print $0}' > ${MATCHES_RAW}

OUTPUT=!{binID}_kmerThreshold_!{params.steps.plasmid.PLSDB.sharedKmerThreshold}.tsv

# extract metadata of matches from plsdb
MATCHES_CONTENT=matches_detailed.tsv
while IFS=$"\t" read line ; do 
	hit=$(echo "$line" | cut -f 1); 
	grep -w $hit ${PLSDB}/plsdb.tsv \
		| sed "s/^/${line}\t/g" \
		| sed "s/^/!{sample}\t!{binID}\t/g" ; 
done <${MATCHES_RAW} > ${MATCHES_CONTENT}

if [ -s ${MATCHES_CONTENT} ]; then
  # add mash query, reference, p-value and distance fields
  cat <(head -n 1 ${PLSDB}/plsdb.tsv \
	| sed "s/^/MASH_REFERENCE\tMASH_QUERY\tMASH_DISTANCE\tMASH_P_VALUE\tMASH_MATCHING_HASHES\t/g"  \
	| sed 's/^/SAMPLE\tBIN_ID\t/g') ${MATCHES_CONTENT} > ${OUTPUT}
fi
