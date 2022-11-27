set -e

OPTIONS=$1
YAML="${2:-example_params/plasmid.yml}"
WORK="${3:-work}"
PROFILE="${4:-standard}"
make run_small_full_test WORK_DIR=${WORK} \
        PARAMS_FILE=$YAML \
       	PROFILE="${PROFILE}" \
       	OPTIONS=" $OPTIONS " \
        ENTRY="wPlasmids"
make check
