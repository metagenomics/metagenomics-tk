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
