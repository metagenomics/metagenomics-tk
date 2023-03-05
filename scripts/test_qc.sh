set -e

ENTRY="wShortReadQualityControl"
OPTIONS=$1
YAML="${2:-example_params/qc.yml}"
WORK="${3:-work}_${ENTRY}"
PROFILE="${4:-standard}"
LOG_DIR="${WORK}/logs"
make run_small_full_test WORK_DIR=${WORK} \
        PARAMS_FILE=$YAML \
       	PROFILE="$PROFILE" \
	LOG_DIR=${LOG_DIR} \
       	OPTIONS=" $OPTIONS " \
        ENTRY="${ENTRY}"

make check LOG_DIR=${LOG_DIR}
