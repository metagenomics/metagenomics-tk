set -e

ENTRY="wReadMapping"
OPTIONS=$1
YAML="${2:-example_params/readMapping.yml}" 
WORK="${3:-work}_${ENTRY}"
PROFILE="${4:-standard}"
LOG_DIR="${WORK}/logs"
make run_small_full_test \
	WORK_DIR=${WORK} OPTIONS=" $OPTIONS " \
       	PROFILE="${PROFILE}"  ENTRY="${ENTRY}" \
	LOG_DIR=${LOG_DIR} \
	PARAMS_FILE=$YAML

make check LOG_DIR=${LOG_DIR}
