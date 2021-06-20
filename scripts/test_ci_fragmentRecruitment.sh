make run_small_full_test CACHE="NXF_HOME=/mnt/nextflow_cache/"  WORK_DIR="/mnt/work/work" OPTIONS=" -process.cache='lenient'   --steps.fragmentRecruitment.frhit.genomes=test/bins/small/bin.*.fa --steps.fragmentRecruitment.frhit.samples=test/reads/small/reads.tsv " PROFILE="local"  ENTRY="wFragmentRecruitment" PARAMS_FILE=example_params/fragmentRecruitment.yml || exit 1
make check
