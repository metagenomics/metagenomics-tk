OPTIONS=$1
YAML="${2:-example_params/readMapping.yml}" 
make run_small_full_test \
	WORK_DIR="work" OPTIONS=" $OPTIONS " \
       	PROFILE="local"  ENTRY="wReadMapping" \
	PARAMS_FILE=$YAML

make check