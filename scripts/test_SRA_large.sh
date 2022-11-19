set -e

bash ./scripts/test_SRA.sh " -c  test_data/assets/aws.config " ./example_params/SRAmany.yml $1

# Check if the produces SRA IDs are exactly the ones that we used as input

set +e

MISSING_DATASETS="$(cat <(cut -f 1 log/samplesILLUMINA.tsv | tail -n +2) <(cat test_data/SRA/manyDatasets.txt | tail -n +2) \
	| sort | uniq -c | grep -v ' 2 ')"

if [ -z "$MISSING_DATASETS" ]
then
	echo "Only the input datasets are reported."
	exit 0
else
	echo "Datasets are missing"
	exit 1
fi
