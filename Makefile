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
	override ENTRY = "wFullPipeline"
endif

ifndef LOG_DIR
	override LOG_DIR = "log"
endif

ifndef OPTIONS
	override OPTIONS = ""
endif

ifndef PARAMS_FILE
	override PARAMS_FILE = ${CURRENT_DIR}/example_params/fullPipeline.yml
endif

ifdef PRESET
	override PARAMS_FILE=""
endif

ifneq (${PARAMS_FILE}, "")
	override PARAMS_COMMAND =  -params-file ${PARAMS_FILE}
endif

ifndef BRANCH
	override BRANCH = "dev"
endif

ifndef SECRET_NAME
	override SECRET_NAME = ""
endif

ifndef SECRET_VALUE
	override SECRET_VALUE = ""
endif

ifndef VERSION
	override VERSION = $(shell cat VERSIONS.txt | sort | tail -n 1)
endif

ifndef MODULE_DB_TEST_EXTRACTED
	override MODULE_DB_TEST_EXTRACTED = "/vol/spool/peter/plsdb"
	override MODULE_DB_TEST_MD5SUM = "caa846fbb689deba0e35ef559793b521"
	override MODULE_DB_TEST_HTTPS = "https://openstack.cebitec.uni-bielefeld.de:8080/databases/plsdb.zip"
	override MODULE_DB_TEST_PATH = "/vol/spool/peter/plsdb.zip"
	override MODULE_DB_TEST_S3PATH = "s3://databases/plsdb.zip"
	override MODULE_DB_TEST_S3_DIRECTORY_PATH = "s3://databases/plsdb/"
	override MODULE_DB_TEST_S5CMD_COMMAND = " --retry-count 30 --no-verify-ssl --endpoint-url https://openstack.cebitec.uni-bielefeld.de:8080 "
	override MODULE_DB_TEST_GENERATED_YML = /vol/spool/peter/meta-omics-toolkit/generated_yamls/
	override MODULE_DB_TEST_YML = example_params/plasmid.yml
	override MODULE_DB_TEST_YML_PATH = ".steps.plasmid.PLSDB.database=env(database)"
	override MODULE_DB_TEST_YML_SCRIPT = "./scripts/test_plasmids.sh" 
	override MODULE_DB_TEST_GENERATED_YML_DIR = "plasmid_yaml_database_tests"
	override MODULE_DB_TEST_SKIP_TESTS = ""
	override MODULE_DB_TEST_REMOVE_DB = "no"
endif

runDatabaseTest: ## Run database tests
	bash ./scripts/test_module_db.sh ${MODULE_DB_TEST_EXTRACTED} ${MODULE_DB_TEST_MD5SUM} \
		${MODULE_DB_TEST_HTTPS} ${MODULE_DB_TEST_PATH} ${MODULE_DB_TEST_S3PATH} ${MODULE_DB_TEST_S3_DIRECTORY_PATH} \
		${MODULE_DB_TEST_S5CMD_COMMAND} ${MODULE_DB_TEST_GENERATED_YML} \
		${MODULE_DB_TEST_YML} ${MODULE_DB_TEST_YML_PATH} ${MODULE_DB_TEST_YML_SCRIPT} ${MODULE_DB_TEST_GENERATED_YML_DIR} \
		${MODULE_DB_TEST_SKIP_TESTS} ${MODULE_DB_TEST_REMOVE_DB}

.PHONY: list clean test_clean run_small_full_test check changelog python_version_check checkPublisDirMode

clean : ## Removes all files that are produced during runs are not necessary
	- rm .nextflow.log*
	- rm report.html*
	- rm timeline.html*
	- rm trace.txt*
	- rm dag.dot*

changelog: ## Creates a new CHANGELOG
	LATEST="$$(git describe --tags $$(git rev-list --tags --max-count=1))"; \
	docker run -v "$${PWD}":/workdir quay.io/git-chglog/git-chglog:latest "$${LATEST}" 

changelogTag: ## Creates a new CHANGELOG by specifying a tag with the environment variable TAG.
	docker run -v "$${PWD}":/workdir quay.io/git-chglog/git-chglog:latest "${TAG}"

nextflow: ## Downloads Nextflow binary
	wget -qO- https://github.com/nextflow-io/nextflow/releases/download/v${VERSION}/nextflow-${VERSION}-all > nextflow
	chmod a+x nextflow

check: ## Checks if processes did fail and the current nextflow run returns exit code 1. (Useful in a github actions context)
	bash ./scripts/check_log.sh $$(ls -1 ${LOG_DIR}/trace.* | tail -n 1 ) || (echo "$?"; exit 1)

checkPublisDirMode: ## Check if publishDirMode is set in process
	! (grep -r publishIR modules/ | grep -v "//" | grep -v params.publishDirMode) || echo "publishDirMode must always be set in process publishDir method"

python_version_check:
	./scripts/pythonVersionCheck.sh

wiki_venv: python_version_check ## Install virtual environment for wiki
	python3 -m venv wiki_venv
	. wiki_venv/bin/activate && pip install -r wiki_scripts/requirements.txt

wiki_publish: wiki_venv ## Publish the local wiki to ghpages
	( \
		. wiki_venv/bin/activate; \
                ./wiki_venv/bin/mike deploy --config-file wiki_scripts/mkdocs.yml --push --update-aliases ${TOOLKIT_TAG} latest; \
	)
	
dev_wiki: wiki_venv ## Start mkdocs developer session
	( \
		. wiki_venv/bin/activate; \
		mkdocs serve --config-file wiki_scripts/mkdocs.yml; \
	)

build_publish_docker:
	bash ./scripts/buildPublishImage.sh ${COMMIT_START} ${COMMIT_END} ${DOCKER_REPOSITORY}
	
set_secrets: nextflow ## Set secrets for sensitive data access
	NXF_HOME=$$PWD/.nextflow ./nextflow secrets set ${SECRET_NAME} ${SECRET_VALUE}

run_small_full_test: nextflow ## Prepares input files like downloading bins and reads and executes Nextflow. The default configuration it runs the full pipeline locally.
	NXF_HOME=$$PWD/.nextflow ./nextflow run main.nf ${OPTIONS} -work-dir ${WORK_DIR} -profile ${PROFILE} -resume -entry ${ENTRY} ${PARAMS_COMMAND} --logDir ${LOG_DIR} ; exit $$?


help: ## Lists available Makefile commands
	@egrep '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | sort | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-16s\033[0m %s\n", $$1, $$2}'
