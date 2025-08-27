set -e

ENTRY="wFullPipeline"
OPTIONS=$1
YAML="${2:-example_params/fullPipeline.yml}"
WORK="${3:-work}_${ENTRY}"
PROFILE="${4:-standard}"
VERSION="${5:-}"
LOG_DIR="${6:-${WORK}/logs}"
make run_small_full_test WORK_DIR=${WORK} \
        PARAMS_FILE=$YAML \
	LOG_DIR=${LOG_DIR} \
       	PROFILE="$PROFILE" \
       	OPTIONS=" $OPTIONS " \
        ENTRY="${ENTRY}" \
        VERSION="$VERSION" 

# In case the log dir is stored in S3, then download first the directory.
if [[ "$LOG_DIR" == s3://* ]]; then
  NEW_LOG_DIR=${WORK}/logs
  mkdir -p ${NEW_LOG_DIR}
  ./bin/s5cmd --endpoint-url https://openstack.cebitec.uni-bielefeld.de:8080 cp ${LOG_DIR}/* ${NEW_LOG_DIR}
  LOG_DIR=${NEW_LOG_DIR}
fi  
  
make check LOG_DIR=${LOG_DIR}

