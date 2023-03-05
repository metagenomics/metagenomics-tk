set -e

WORK_DIR=$1

PROFILE=$2

ENTRY="wSRATable"

bash ./scripts/test_SRA.sh " -c  test_data/assets/aws.config " ./example_params/SRAmany.yml ${WORK_DIR} ${PROFILE} 

# Check if the produces SRA IDs are exactly the ones that we used as input

set +e

MISSING_DATASETS="$(cat <(cut -f 1 "${WORK_DIR}_${ENTRY}/logs/samplesILLUMINA.tsv" | tail -n +2) <(cat test_data/SRA/manyDatasets.txt | tail -n +2) \
	| sort | uniq -c | grep -v ' 2 ')"

if [ -z "$MISSING_DATASETS" ]
then
	echo "Only the input datasets are reported."
	exit 0
else
	echo "Datasets are missing"
	exit 1
fi
