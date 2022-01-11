set -e

OPTIONS=$1
YAML="${2:-example_params/SRA.yml}" 
PROFILE="${3:-standard}" 
make run_small_full_test \
	WORK_DIR="work" OPTIONS=" $OPTIONS " \
       	PROFILE="$PROFILE"  ENTRY="wSRATable" \
	PARAMS_FILE=$YAML
make check
