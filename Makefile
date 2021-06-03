CURRENT_DIR = $(shell pwd)

ifndef PROFILE
	override PROFILE = "local,conda"
endif

ifndef WORK_DIR
	override WORK_DIR = "work"
endif



small_reads_folder = ${CURRENT_DIR}/test/reads/small
reads_split_test = ${small_reads_folder}/reads_split.tsv
small_read1 = ${small_reads_folder}/read1_1.fq.gz 
small_read2 = ${small_reads_folder}/read1_1.fq.gz 
full_run = ${CURRENT_DIR}/example_params/full_pipeline_params.yml

.PHONY : clean init_test test_clean run_small_full_test
clean :
	- rm .nextflow.log*
	- rm report.html*
	- rm timeline.html*
	- rm trace.txt*
	- rm dag.dot*

nextflow:
	- wget -qO- https://get.nextflow.io | bash

test_clean:
	- rm -rf test/reads

test/reads: test_clean
	- mkdir -p test/reads/small
	- wget https://openstack.cebitec.uni-bielefeld.de:8080/swift/v1/meta_test/small/read1_1.fq.gz -P test/reads/small/
	- wget https://openstack.cebitec.uni-bielefeld.de:8080/swift/v1/meta_test/small/read2_1.fq.gz -P test/reads/small/
	- echo "SAMPLE\tREADS1\tREADS2" > ${reads_split_test}
	- echo "test1\t${small_read1}\t${small_read2}" >> ${reads_split_test}
	- echo "test2\t${small_read1}\t${small_read2}" >> ${reads_split_test}

run_small_full_test: init_test
	- ./nextflow run main.nf -work-dir ${WORK_DIR} -profile ${PROFILE} -resume -entry wPipeline -params-file ${full_run} 
