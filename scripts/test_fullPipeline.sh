OPTIONS=$@
YAML="${2:-example_params/fullPipeline.yml}"
make run_small_full_test WORK_DIR="work" \
        PARAMS_FILE=$YAML \
       	PROFILE="local" \
       	OPTIONS=" $OPTIONS " \
        ENTRY="wPipeline"

make check
