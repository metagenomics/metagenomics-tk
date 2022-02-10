CURRENT_DIR = $(shell pwd)


ifndef PROFILE
	override PROFILE = "standard,conda"
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

ifndef BRANCH
	override BRANCH = "dev"
endif


ifndef MODULE_DB_TEST_EXTRACTED
	override MODULE_DB_TEST_EXTRACTED = "/vol/spool/peter/plsdb"
	override MODULE_DB_TEST_MD5SUM = "caa846fbb689deba0e35ef559793b521"
	override MODULE_DB_TEST_HTTPS = "https://openstack.cebitec.uni-bielefeld.de:8080/databases/plsdb.zip"
	override MODULE_DB_TEST_PATH = "/vol/spool/peter/plsdb.zip"
	override MODULE_DB_TEST_S3PATH = "s3://databases/plsdb.zip"
	override MODULE_DB_TEST_S3_DIRECTORY_PATH = "s3://databases/plsdb/"
	override MODULE_DB_TEST_S5CMD_COMMAND = " --retry-count 30 --no-verify-ssl --endpoint-url https://openstack.cebitec.uni-bielefeld.de:8080 "
	override MODULE_DB_TEST_CREDENTIALS = /vol/spool/peter/meta-omics-toolkit/credentials
	override MODULE_DB_TEST_GENERATED_YML = /vol/spool/peter/meta-omics-toolkit/generated_yamls/
	override MODULE_DB_TEST_YML = example_params/plasmid.yml
	override MODULE_DB_TEST_YML_PATH = ".steps.plasmid.PLSDB.database=env(database)"
	override MODULE_DB_TEST_YML_SCRIPT = "./scripts/test_plasmids.sh" 
	override MODULE_DB_TEST_GENERATED_YML_DIR = "plasmid_yaml_database_tests"
	override MODULE_DB_TEST_SKIP_TESTS = ""
endif

runDatabaseTest: ## Run database tests
	bash ./scripts/test_module_db.sh ${MODULE_DB_TEST_EXTRACTED} ${MODULE_DB_TEST_MD5SUM} \
		${MODULE_DB_TEST_HTTPS} ${MODULE_DB_TEST_PATH} ${MODULE_DB_TEST_S3PATH} ${MODULE_DB_TEST_S3_DIRECTORY_PATH} \
		${MODULE_DB_TEST_S5CMD_COMMAND} ${MODULE_DB_TEST_CREDENTIALS} ${MODULE_DB_TEST_GENERATED_YML} \
		${MODULE_DB_TEST_YML} ${MODULE_DB_TEST_YML_PATH} ${MODULE_DB_TEST_YML_SCRIPT} ${MODULE_DB_TEST_GENERATED_YML_DIR} \
		${MODULE_DB_TEST_SKIP_TESTS}

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
	! grep -q "FAILED" log/trace.tsv || (echo "$?"; exit 1)


wiki_venv: ## Install virtual environment for wiki
	python3 -m venv wiki_venv
	. wiki_venv/bin/activate && pip install -r wiki_scripts/requirements.txt

wiki_build.html: wiki_venv ## Build wiki html
	( \
		. wiki_venv/bin/activate; \
		cd wiki_scripts; \
		mkdocs build; \
		htmlark site/print_page/index.html -o ../wiki_build.html \
	)

publish_wiki_html: wiki_build.html ## Publish wiki html
	mc cp wiki_build.html dereplication/meta-omics-toolkit/${BRANCH}.html
	
dev_wiki: wiki_venv ## Start mkdocs developer session
	( \
		. wiki_venv/bin/activate; \
		cd wiki_scripts; \
		mkdocs serve; \
	)

	
run_small_full_test: nextflow ## Prepares input files like downloading bins and reads and executes Nextflow. The default configuration it runs the full pipeline locally.
	./nextflow run main.nf ${OPTIONS} -work-dir ${WORK_DIR}_${ENTRY} -profile ${PROFILE} -resume -entry ${ENTRY} -params-file ${PARAMS_FILE}; exit $$?


help: ## Lists available Makefile commands
	@egrep '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | sort | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-16s\033[0m %s\n", $$1, $$2}'
