make run_small_full_test WORK_DIR="/mnt/work/work" OPTIONS=" --steps.magAttributes.input=test/bins/small/attributes.tsv " PROFILE="local"  ENTRY="wMagAttributes" PARAMS_FILE=example_params/magAttributes.yml  || exit 1
make check
