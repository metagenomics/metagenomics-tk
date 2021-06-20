make run_small_full_test DEST="/mnt/nextflow_cache" \
       	WORK_DIR="/mnt/work/work"  PROFILE="local" \
       	OPTIONS=" --inputs=/mnt/nextflow_cache/test/reads/small/reads_split.tsv --steps.magAttributes.checkm.database=/mnt/checkm --steps.magAttributes.gtdb.database=/mnt/gtdb/release202 " || exit 1
make check
