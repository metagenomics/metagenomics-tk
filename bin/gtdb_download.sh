EXTRACTED_DB=$1
DOWNLOAD_LINK=$2
S5CMD_PARAMS=$3
CPUS=$4
POLISHED_DB=$5
MD5SUM=$6


# Check developer documentation
GTDB=""
if [ -z "${EXTRACTED_DB}" ]
then
     DATABASE=${POLISHED_DB}gtdb
     LOCK_FILE=${DATABASE}/lock.txt

     # Download gtdb if necessary
     mkdir -p ${DATABASE}
     flock ${LOCK_FILE} concurrentDownload.sh --output=${DATABASE} \
   	--link=${DOWNLOAD_LINK} \
   	--httpsCommand="wget -O gtdb.tar.gz ${DOWNLOAD_LINK} && tar xzvf gtdb.tar.gz && rm gtdb.tar.gz" \
   	--s3FileCommand="s5cmd ${S5CMD_PARAMS} cp --concurrency ${CPUS} ${DOWNLOAD_LINK} gtdb.tar.gz  && tar xzvf gtdb.tar.gz && rm gtdb.tar.gz " \
           --s3DirectoryCommand="s5cmd ${S5CMD_PARAMS} cp --concurrency ${CPUS} ${DOWNLOAD_LINK} . " \
   	--s5cmdAdditionalParams="${S5CMD_PARAMS}" \
   	--localCommand="tar -xzvf ${DOWNLOAD_LINK} " \
   	--expectedMD5SUM=${MD5SUM}

     GTDB=$(readlink -f ${DATABASE}/out/*)
else
     GTDB=${EXTRACTED_DB}
fi

echo $GTDB > gtdbPath.txt