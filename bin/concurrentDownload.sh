#!/bin/bash

while [ $# -gt 0 ]; do
	  case "$1" in
	    --output=*) OUTPUT_PATH="${1#*=}"
	    ;;
	    --checkpoint=*) CHECKPOINT_FILE="${1#*=}"
	    ;;
	    --expectedVersion=*) EXPECTED="${1#*=}"
	    ;;
	    --expectedMD5SUM=*) EXPECTED="${1#*=}"
	    ;;
	    --command=*) COMMAND="${1#*=}"
	    ;;
	    --mode=*) MODE="${1#*=}"
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
mkdir -p ${DATABASE_OUT}

compare(){
   MODE=$1
   case $MODE in
	MD5SUM)
      echo "MD5SUM"
      cd ${DATABASE_OUT}
      ls -1
      MD5SUM=$(find . -type f -exec md5sum {} + | sort | md5sum | cut -d ' ' -f 1)
      T=$(test "${MD5SUM}" = "${EXPECTED}")
      echo "${MD5SUM} md5SUM does not match, ${T}"
      test "${MD5SUM}" = "${EXPECTED}" && grep -Fxq ${EXPECTED} $CHECKPOINT_FILE
      ;;
        VERSION)
      echo "VERSION";
      grep -Fxq $EXPECTED $CHECKPOINT_FILE;
      ;;
      *)
      echo "Please provide either MD5SUM or VERSION as mode."
      ;;
   esac
   return $?
}

if [ -f "$CHECKPOINT_FILE" ] && compare $MODE ; then
    echo "Database already exists!"
else 
    echo "Database will be downloaded!"
    cd ${DATABASE_OUT}
    eval "$COMMAND"
    echo ${EXPECTED} > ${CHECKPOINT_FILE}
fi
