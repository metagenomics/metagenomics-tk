OPTIONS=$@
make run_small_full_test \
	WORK_DIR="work" OPTIONS=" $OPTIONS " \
       	PROFILE="local"  ENTRY="wMagAttributes" \
	PARAMS_FILE=example_params/magAttributes.yml

make check
