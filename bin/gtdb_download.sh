EXTRACTED_DB=$1
DOWNLOAD_LINK=$2
S5CMD_PARAMS=$3
CPUS=$4
POLISHED_DB=$5
MD5SUM=$6
S3_gtdb_ACCESS=$7
S3_gtdb_SECRET=$8


# Check developer documentation
GTDB=""
if [ -z "${EXTRACTED_DB}" ]
then
     DATABASE=${POLISHED_DB}gtdb
     LOCK_FILE=${DATABASE}/lock.txt

     if [ ! -z "${S3_gtdb_ACCESS}" ]
     then
          export AWS_ACCESS_KEY_ID=${S3_gtdb_ACCESS}
          export AWS_SECRET_ACCESS_KEY=${S3_gtdb_SECRET}
     fi

     # Download gtdb if necessary
     mkdir -p ${DATABASE}
     flock ${LOCK_FILE} concurrentDownload.sh --output=${DATABASE} \
   	--link=${DOWNLOAD_LINK} \
   	--httpsCommand="wgetStatic -qO- ${DOWNLOAD_LINK} | tar xvz " \
   	--s3FileCommand="s5cmd ${S5CMD_PARAMS} cat ${DOWNLOAD_LINK} | tar xzv " \
           --s3DirectoryCommand="s5cmd ${S5CMD_PARAMS} cp ${DOWNLOAD_LINK} . " \
   	--s5cmdAdditionalParams="${S5CMD_PARAMS}" \
   	--localCommand="tar -xzvf ${DOWNLOAD_LINK} " \
   	--expectedMD5SUM=${MD5SUM}

     GTDB=$(readlink -f ${DATABASE}/out/*)
else
     GTDB=${EXTRACTED_DB}
fi

echo $GTDB > gtdbPath.txt
