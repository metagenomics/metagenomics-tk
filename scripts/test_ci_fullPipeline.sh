make run_small_full_test CACHE="NXF_HOME=/mnt/nextflow_cache/" WORK_DIR="/mnt/work/work"  PROFILE="local" OPTIONS=" --steps.magAttributes.checkm.database=/mnt/checkm -process.cache='lenient'  --steps.magAttributes.gtdb.database=/mnt/gtdb/release202 " || exit 1
make check
