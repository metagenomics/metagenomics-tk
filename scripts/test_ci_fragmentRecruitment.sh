make run_small_full_test DEST="/mnt/nextflow_cache" \
       	WORK_DIR="/mnt/work/work" \
       	OPTIONS=" --steps.fragmentRecruitment.frhit.genomes=/mnt/nextflow_cache/test/bins/small/bin.*.fa --steps.fragmentRecruitment.frhit.samples=/mnt/nextflow_cache/test/reads/small/reads.tsv " \
       	PROFILE="local"  ENTRY="wFragmentRecruitment" PARAMS_FILE=example_params/fragmentRecruitment.yml || exit 1
make check
