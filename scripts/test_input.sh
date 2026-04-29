set -e

PARAMETERS=$1

EXPECTED_SAMPLES=$2

WORK_DIR=$3

PROFILE=$4


ENTRY="wOutputTable"

# Cleanup existing output tables
rm ${WORK_DIR}_${ENTRY}/logs/samples*.tsv

bash ./scripts/test_SRA.sh "${PARAMETERS}" "./example_params/input.yml" ${WORK_DIR} ${PROFILE} 

# Check if the produces SRA IDs are exactly the ones that we used as input

set +e

DIFF="$(diff <(sort -k1,1 ${WORK_DIR}_${ENTRY}/logs/samples*.tsv) <(sort -k 1,1 ${EXPECTED_SAMPLES}))"

if [ -z "$DIFF" ]
then
	echo "Only the input datasets are reported."
	exit 0
else
	echo "The output was not expected."
	echo ${DIFF}
	exit 1
fi
