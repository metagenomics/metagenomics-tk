set -e

ENTRY="wSaveSettings"
OPTIONS=$1
YAML="${2:-example_params/fullPipeline.yml}"
WORK="${3:-work}_${ENTRY}"
PROFILE="${4:-standard}"
VERSION="${5:-}"
LOG_DIR="${WORK}/logs"
make run_small_full_test WORK_DIR=${WORK} \
	        PARAMS_FILE=$YAML \
		LOG_DIR=${LOG_DIR} \
		PROFILE="$PROFILE" \
	        OPTIONS=" $OPTIONS " \
	        ENTRY="${ENTRY}" \
	        VERSION="$VERSION" 

make check LOG_DIR=${LOG_DIR}
