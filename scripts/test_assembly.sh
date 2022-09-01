set -e

OPTIONS=$1
YAML="${2:-example_params/assembly.yml}" 
WORK="${3:-work}"
PROFILE="${4:-standard}"
PID_PATH="${5:-}"
make run_small_full_test \
	WORK_DIR=${WORK} OPTIONS=" $OPTIONS " \
        PID_PATH="$PID_PATH" \
       	PROFILE="$PROFILE"  ENTRY="wAssembly" \
	PARAMS_FILE=$YAML
make check
