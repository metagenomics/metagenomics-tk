OPTIONS=$@
make run_small_full_test WORK_DIR="work" OPTIONS=" $OPTIONS " PROFILE="local"  ENTRY="wDereplication" PARAMS_FILE=example_params/dereplication.yml
make check
