OPTIONS=$@
make run_small_full_test WORK_DIR="work" PARAMS_FILE="example_params/fullPipelineQC.yml"  PROFILE="local" OPTIONS=" $OPTIONS "
make check
