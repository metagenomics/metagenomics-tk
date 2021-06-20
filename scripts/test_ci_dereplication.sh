make run_small_full_test CACHE="NXF_HOME=/mnt/nextflow_cache/"  WORK_DIR="/mnt/work/work" OPTIONS=" -process.cache='lenient' --steps.dereplication.pasolli.input=test/bins/small/attributes.tsv " PROFILE="local"  ENTRY="wDereplication" PARAMS_FILE=example_params/dereplication.yml || exit 1
make check
