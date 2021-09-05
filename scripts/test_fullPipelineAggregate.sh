OPTIONS=$@
YAML="example_params/fullPipelineAggregate.yml"
make run_small_full_test WORK_DIR="work" \
	  PARAMS_FILE=$YAML \
	  PROFILE="local" \
	  OPTIONS=" $OPTIONS " \
	  ENTRY="wAggregatePipeline"
make check
