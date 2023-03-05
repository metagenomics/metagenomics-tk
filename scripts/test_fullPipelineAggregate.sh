set -e

ENTRY="wAggregatePipeline"
OPTIONS=$1
YAML="${2:-example_params/fullPipelineAggregate.yml}"
WORK="${3:-work}_${ENTRY}"
PROFILE="${4:-standard}"
LOG_DIR="${WORK}/logs"
make run_small_full_test WORK_DIR=${WORK} \
	  PARAMS_FILE=$YAML \
	  PROFILE=" $PROFILE " \
	  LOG_DIR=${LOG_DIR} \
	  OPTIONS=" $OPTIONS " \
	  ENTRY="${ENTRY}"
make check LOG_DIR=${LOG_DIR}
