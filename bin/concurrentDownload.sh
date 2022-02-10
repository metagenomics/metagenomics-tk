#!/bin/bash

while [ $# -gt 0 ]; do
	  case "$1" in
	    --output=*) OUTPUT_PATH="${1#*=}"
	    ;;
	    --expectedMD5SUM=*) EXPECTED_MD5SUM="${1#*=}"
	    ;;
	    --s3FileCommand=*) S3_FILE_COMMAND="${1#*=}"
	    ;;
	    --s3DirectoryCommand=*) S3_DIRECTORY_COMMAND="${1#*=}"
	    ;;
	    --s5cmdAdditionalParams=*) S5CMD_ADDITIONAL_PARAMS="${1#*=}"
	    ;;
	    --httpsCommand=*) HTTPS_COMMAND="${1#*=}"
	    ;;
	    --localCommand=*) LOCAL_COMMAND="${1#*=}"
	    ;;
	    --link=*) LINK="${1#*=}"
	    ;;
	    *)
              printf "***************************\n"
              printf "* Error: Invalid argument.*\n"
	      printf "***************************\n"
	      exit 1
	    ;;
	  esac
          shift
done

DATABASE_OUT=${OUTPUT_PATH}/out
MD5SUM_FILE=${OUTPUT_PATH}/md5sum.txt
mkdir -p ${DATABASE_OUT}

function getCommand() {
    if [[ $LINK == s3://* ]]
    then
	if [[ $(s5cmd ${S5CMD_ADDITIONAL_PARAMS}  ls ${LINK} | wc -l) == 1 ]]; then
		echo "$S3_FILE_COMMAND"
	else
		echo "$S3_DIRECTORY_COMMAND"
	fi
    elif [[ $LINK == https://* ]]
    then
    	echo "$HTTPS_COMMAND";
    elif [[ $LINK == /* ]]
    then
    	echo "$LOCAL_COMMAND";
    fi
}

# Compares the expected MD5SUM to the one saved in checkpoint file.
function compareExpectedToCheckpoint() {
   if grep -Fxq ${EXPECTED_MD5SUM} $MD5SUM_FILE; then
     return 0
   else
     return 1
   fi
}

# Retrieve MD5SUM based on MD5SUM of MD5SUMs
function getMD5SUM() {
   cd ${DATABASE_OUT}
   MD5SUM=$(find * -type f -exec md5sum {} + | sort | cut -d ' ' -f 1 | md5sum | cut -d ' ' -f 1)
   echo ${MD5SUM}
}

if [ -f "$MD5SUM_FILE" ] && compareExpectedToCheckpoint $MODE ; then
    echo "Database already exists!"
else 
    echo "Database will be downloaded!"
    # remove wrong version of the database
    rm -rf ${DATABASE_OUT} ${MD5SUM_FILE}

    mkdir -p ${DATABASE_OUT}

    # download database
    cd ${DATABASE_OUT}
    COMMAND=$(getCommand)
    eval "$COMMAND"        

    # compute MD5SUM and save it in the checkpoint file
    MD5SUM=$(getMD5SUM)
    echo "MD5SUM: ${MD5SUM}"
    echo "EXPECTED MD5SUM: ${EXPECTED_MD5SUM}"
    if [[ "${MD5SUM}" == "${EXPECTED_MD5SUM}" ]]; then
            echo ${MD5SUM} > ${MD5SUM_FILE}
    else
	    echo " Computed MD5SUM does not match to the expected one. "
	    exit 1
    fi
fi
