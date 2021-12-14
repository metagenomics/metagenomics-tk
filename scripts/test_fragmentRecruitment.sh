set -e

OPTIONS=$1
YAML="${2:-example_params/fragmentRecruitment.yml}"
WORK="${3:-work}"
make run_small_full_test WORK_DIR=${WORK} \
	OPTIONS=" $OPTIONS  " \
	PROFILE="standard" \
	ENTRY="wFragmentRecruitment" \
	PARAMS_FILE=${YAML}
make check
