set -e

OPTIONS=$1
YAML="${2:-example_params/fullPipelineAggregate.yml}"
WORK="${3:-work}"
PROFILE="${4:-standard}"
PID_PATH="${5:-}"
make run_small_full_test WORK_DIR=${WORK} \
	  PARAMS_FILE=$YAML \
          PID_PATH="$PID_PATH" \
	  PROFILE=" $PROFILE " \
	  OPTIONS=" $OPTIONS " \
	  ENTRY="wAggregatePipeline"
make check
