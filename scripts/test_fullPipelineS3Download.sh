OPTIONS=$1
YAML="${2:-example_params/fullPipelineQC.yml}"
WORK="${3:-work}"
make run_small_full_test \ 
	WORK_DIR=${WORK} \
	PARAMS_FILE=${YAML} \
	PROFILE="local" \
	OPTIONS=" $OPTIONS "
make check
