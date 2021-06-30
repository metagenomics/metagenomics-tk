make run_small_full_test DEST="/mnt/nextflow_cache" \
       	WORK_DIR="/mnt/work/work"  PROFILE="local" \
       	OPTIONS=" --input=/mnt/nextflow_cache/test/reads/small/reads_split.tsv --steps.magAttributes.checkm.database=/mnt/checkm --steps.magAttributes.gtdb.database=/mnt/gtdb/release202 " \
        PARAMS_FILE="example_params/fullPipeline_fraction/fullPipeline_fraction_magAttributes.yml" || exit 1
make check
