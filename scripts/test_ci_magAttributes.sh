make run_small_full_test CACHE="NXF_HOME=/mnt/nextflow_cache/" WORK_DIR="/mnt/work/work" OPTIONS=" -process.cache='lenient' --steps.magAttributes.input=test/bins/small/attributes.tsv " PROFILE="local"  ENTRY="wMagAttributes" PARAMS_FILE=example_params/magAttributes.yml  || exit 1
make check
