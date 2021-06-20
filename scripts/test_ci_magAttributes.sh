make run_small_full_test DEST="/mnt/nextflow_cache" WORK_DIR="/mnt/work/work" OPTIONS="  --steps.magAttributes.checkm.database=/mnt/checkm  --steps.magAttributes.gtdb.database=/mnt/gtdb/release202  --steps.magAttributes.input=test/bins/small/attributes.tsv " PROFILE="local"  ENTRY="wMagAttributes" PARAMS_FILE=example_params/magAttributes.yml  || exit 1
make check
