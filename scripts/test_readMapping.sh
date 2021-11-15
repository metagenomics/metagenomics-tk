set -e

OPTIONS=$1
YAML="${2:-example_params/readMapping.yml}" 
WORK="${3:-work}"
make run_small_full_test \
	WORK_DIR=${WORK} OPTIONS=" $OPTIONS " \
       	PROFILE="local"  ENTRY="wReadMapping" \
	PARAMS_FILE=$YAML

make check
