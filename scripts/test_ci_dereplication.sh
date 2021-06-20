make run_small_full_test DEST="/mnt/nextflow_cache"  WORK_DIR="/mnt/work/work" OPTIONS="  --steps.dereplication.pasolli.input=test/bins/small/attributes.tsv " PROFILE="local"  ENTRY="wDereplication" PARAMS_FILE=example_params/dereplication.yml || exit 1
make check
