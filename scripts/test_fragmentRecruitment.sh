set -e

ENTRY="wFragmentRecruitment"
OPTIONS=$1
YAML="${2:-example_params/fragmentRecruitment.yml}"
WORK="${3:-work}_${ENTRY}"
PROFILE="${4:-standard}"
LOG_DIR="${WORK}/logs"
make run_small_full_test WORK_DIR=${WORK} \
	OPTIONS=" $OPTIONS  " \
	PROFILE="$PROFILE" \
	LOG_DIR=${LOG_DIR} \
	ENTRY="${ENTRY}" \
	PARAMS_FILE=${YAML}
make check LOG_DIR=${LOG_DIR}
