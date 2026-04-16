set -e

WORK_DIR=$1

PROFILE=$2

ENTRY="wOutputTable"

bash ./scripts/test_SRA.sh " -c  test_data/assets/aws.config \
				      --input.SRA.NCBI.id 'SRR5651439' \
                                      --input.SRA.S3.id 'ERR12263778 SRR29912082' \
                                      --input.ont.r https://openstack.cebitec.uni-bielefeld.de:8080/swift/v1/meta_test/SRR16328449_qc.fq.gz \
                                      --input.ont.names 'nano' \
                                      --input.paired.r1 'https://openstack.cebitec.uni-bielefeld.de:8080/swift/v1/meta_test/small/read1_1bla.fq.gz https://openstack.cebitec.uni-bielefeld.de:8080/swift/v1/meta_test/small/read1_1.fq.gz' \
                                      --input.paired.r2 'https://openstack.cebitec.uni-bielefeld.de:8080/swift/v1/meta_test/small/read2_1bla.fq.gz https://openstack.cebitec.uni-bielefeld.de:8080/swift/v1/meta_test/small/read2_1.fq.gz' \
                                      --input.paired.names 'test1bla test2' " \
                                      "./example_params/input.yml" ${WORK_DIR} ${PROFILE} 

# Check if the produces SRA IDs are exactly the ones that we used as input

set +e

cat <(cut -f 1 work_wOutputTable/logs/samples*.tsv | grep -v SAMPLE)

MISSING_DATASETS="$(cat <(cut -f 1 ${WORK_DIR}_${ENTRY}/logs/samples*.tsv | grep -v SAMPLE) <(cat test_data/input/expectedSamples.tsv) \
	| sort | uniq -c | grep -v ' 2 ')"

if [ -z "$MISSING_DATASETS" ]
then
	echo "Only the input datasets are reported."
	exit 0
else
	echo "Datasets are missing"
	exit 1
fi
