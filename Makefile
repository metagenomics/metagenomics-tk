CURRENT_DIR = $(shell pwd)


ifndef PROFILE
	override PROFILE = "local,conda"
endif

ifndef DEST
	override DEST = $(shell pwd)
endif

ifndef WORK_DIR
	override WORK_DIR = "work"
endif

ifndef ENTRY
	override ENTRY = "wPipeline"
endif

ifndef OPTIONS
	override OPTIONS = ""
endif

ifndef PARAMS_FILE
	override PARAMS_FILE = ${CURRENT_DIR}/example_params/fullPipeline.yml
endif

small_reads_folder = ${DEST}/test/reads/small
bins_folder = ${DEST}/test/bins/small
reads_split_test = ${small_reads_folder}/reads_split.tsv
reads_interleaved_test = ${small_reads_folder}/reads.tsv
bins_attributes_test = ${bins_folder}/attributes.tsv
small_read1 = ${small_reads_folder}/read1_1.fq.gz 
small_read2 = ${small_reads_folder}/read2_1.fq.gz 
small_read_interleaved = ${small_reads_folder}/interleaved.fq.gz
bin1 = ${bins_folder}/bin.1.fa
bin2 = ${bins_folder}/bin.2.fa
bin3 = ${bins_folder}/bin.8.fasta
bin4 = ${bins_folder}/bin.9.fasta
bin5 = ${bins_folder}/bin.32.fa



.PHONY: list clean test_clean run_small_full_test check changelog
clean : ## Removes all files that are produced during runs are not necessary
	- rm .nextflow.log*
	- rm report.html*
	- rm timeline.html*
	- rm trace.txt*
	- rm dag.dot*

changelog: ## Creates a new CHANGELOG.md file
	LATEST="$$(git describe --tags $$(git rev-list --tags --max-count=1))"; \
	docker run -v "$${PWD}":/workdir quay.io/git-chglog/git-chglog:latest "$${LATEST}"

nextflow: ## Downloads Nextflow binary
	- wget -qO- https://get.nextflow.io | bash

check: ## Checks if processes did failed in the current nextflow returns exit code 1. (Useful in github actions context)
	! grep -q "FAILED" log/trace.tsv && grep -q "Succeeded" nextflow.stdout.log || (echo "$?"; exit 1)

${DEST}/test/reads/small: ## Downloads split fastq files and creates tsv files which can be used as input for meta-omics-toolkit
	- mkdir -p  ${DEST}/test/reads/small
	- wget -q https://openstack.cebitec.uni-bielefeld.de:8080/swift/v1/meta_test/small/read1_1.fq.gz -P ${DEST}/test/reads/small/
	- wget -q https://openstack.cebitec.uni-bielefeld.de:8080/swift/v1/meta_test/small/read2_1.fq.gz -P ${DEST}/test/reads/small/
	- echo "SAMPLE\tREADS1\tREADS2" > ${reads_split_test}
	- echo "test1\t${small_read1}\t${small_read2}" >> ${reads_split_test}
	- echo "test2\t${small_read1}\t${small_read2}" >> ${reads_split_test}

${DEST}/test/reads/small/interleaved.fq.gz: ## Downloads interleaved fastq files and creates tsv files which can be used as input for meta-omics-toolkit
	- mkdir -p ${DEST}/test/reads/small
	- wget -O ${DEST}/test/reads/small/interleaved.fq.gz  -q https://openstack.cebitec.uni-bielefeld.de:8080/swift/v1/meta_test/small/RL_S001__insert_270_new3.fq.gz
	- echo "SAMPLE\tREADS" > ${reads_interleaved_test}
	- echo "test1\t${small_read_interleaved}" >> ${reads_interleaved_test}
	- echo "test2\t${small_read_interleaved}" >> ${reads_interleaved_test}

${DEST}/test/bins/small/: ## Downloads bins and creates a tsv file with bin properties like contamination, completeness.. etc.
	- mkdir -p ${DEST}/test/bins/small
	- wget -O ${DEST}/test/bins/small/bin.1.fa -q https://openstack.cebitec.uni-bielefeld.de:8080/swift/v1/meta_test/bins/bin.1.fa
	- wget -O ${DEST}/test/bins/small/bin.2.fa -q https://openstack.cebitec.uni-bielefeld.de:8080/swift/v1/meta_test/bins/bin.2.fa
	- wget -O ${DEST}/test/bins/small/bin.32.fa -q https://openstack.cebitec.uni-bielefeld.de:8080/swift/v1/meta_test/bins/bin.32.fa
	- wget -O ${DEST}/test/bins/small/bin.8.fasta -q https://openstack.cebitec.uni-bielefeld.de:8080/swift/v1/meta_test/bins/bin.8.fasta
	- wget -O ${DEST}/test/bins/small/bin.9.fasta -q https://openstack.cebitec.uni-bielefeld.de:8080/swift/v1/meta_test/bins/bin.9.fasta
	- echo "DATASET\tBIN_ID\tPATH\tCOMPLETENESS\tCONTAMINATION\tCOVERAGE\tN50\tHETEROGENEITY" >> ${bins_attributes_test}
	- echo "test1\tbin.1\t${bin1}\t100\t0\t10\t5000\t10" >> ${bins_attributes_test}
	- echo "test1\tbin.2\t${bin2}\t100\t0\t10\t5000\t10" >> ${bins_attributes_test}
	- echo "test1\tbin.8\t${bin3}\t100\t0\t10\t5000\t10" >> ${bins_attributes_test}
	- echo "test2\tbin.9\t${bin4}\t100\t0\t10\t5000\t10" >> ${bins_attributes_test}
	- echo "test2\tbin.32\t${bin5}\t100\t0\t10\t5000\t10" >> ${bins_attributes_test}


run_small_full_test: ${DEST}/test/reads/small nextflow ${DEST}/test/bins/small/ ${DEST}/test/reads/small/interleaved.fq.gz ## Prepares input files like downloading bins and reads and executes Nextflow. The default configuration it runs the full pipeline locally.
	./nextflow run main.nf ${OPTIONS} -work-dir ${WORK_DIR}_${ENTRY} -profile ${PROFILE} -resume -entry ${ENTRY} -params-file ${PARAMS_FILE} | tee nextflow.stdout.log 


help: ## Lists available Makefile commands
	@egrep '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | sort | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-16s\033[0m %s\n", $$1, $$2}'
