set -e

OPTIONS=$1
YAML="${2:-example_params/fullPipeline.yml}"
WORK="${3:-work}"
PROFILE="${4:-standard}"
VERSION="${5:-}"
make run_small_full_test WORK_DIR=${WORK} \
        PARAMS_FILE=$YAML \
       	PROFILE="$PROFILE" \
       	OPTIONS=" $OPTIONS " \
        ENTRY="wFullPipeline"
        VERSION="$VERSION" \
make check
