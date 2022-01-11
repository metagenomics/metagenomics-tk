set -e

OPTIONS=$1
YAML="${2:-example_params/dereplication.yml}" 
WORK="${3:-work}"
PROFILE="${4:-standard}"
make run_small_full_test \
	WORK_DIR=${WORK} OPTIONS=" $OPTIONS " \
       	PROFILE="$PROFILE"  ENTRY="wDereplication" \
	PARAMS_FILE=$YAML

make check
