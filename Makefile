CURRENT_DIR = $(shell pwd)


small_reads_folder = ${CURRENT_DIR}/test/reads/small
reads_split_test = ${small_reads_folder}/reads_split.tsv
small_read1 = ${small_reads_folder}/read1_1.fq.gz 
small_read2 = ${small_reads_folder}/read2_1.fq.gz 

ifndef PROFILE
	override PROFILE = "local,conda"
endif

ifndef WORK_DIR
	override WORK_DIR = "work"
endif

ifndef ENTRY
	override ENTRY = "wPipeline"
endif

ifndef PARAMS_FILE
	override PARAMS_FILE = ${CURRENT_DIR}/example_params/full_pipeline_params.yml
endif

.PHONY: list clean test_clean run_small_full_test
clean :
	- rm .nextflow.log*
	- rm report.html*
	- rm timeline.html*
	- rm trace.txt*
	- rm dag.dot*

nextflow:
	- wget -qO- https://get.nextflow.io | bash

check:
	- ! grep -q "FAILED" log/trace.tsv


test/reads/small:
	- mkdir -p test/reads/small
	- wget -q https://openstack.cebitec.uni-bielefeld.de:8080/swift/v1/meta_test/small/read1_1.fq.gz -P test/reads/small/
	- wget -q https://openstack.cebitec.uni-bielefeld.de:8080/swift/v1/meta_test/small/read2_1.fq.gz -P test/reads/small/
	- echo "SAMPLE\tREADS1\tREADS2" > ${reads_split_test}
	- echo "test1\t${small_read1}\t${small_read2}" >> ${reads_split_test}
	- echo "test2\t${small_read1}\t${small_read2}" >> ${reads_split_test}

run_small_full_test: test/reads/small nextflow
	- ./nextflow run main.nf -work-dir ${WORK_DIR}_${ENTRY} -profile ${PROFILE} -resume -entry ${ENTRY} -params-file ${PARAMS_FILE} 

list:
	@$(MAKE) -pRrq -f $(lastword $(MAKEFILE_LIST)) : 2>/dev/null | awk -v RS= -F: '/^# File/,/^# Finished Make data base/ {if ($$1 !~ "^[#.]") {print $$1}}' | sort | egrep -v -e '^[^[:alnum:]]' -e '^$@$$'

