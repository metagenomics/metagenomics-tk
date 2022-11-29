
# create input, output files and run default gtdbtk command
mkdir output
ls -1 !{bins} | xargs -I {} readlink -f {} > bin.path
paste -d$'\t' bin.path <(for p in $(cat bin.path); do basename $p; done) > input.tsv

# Check developer documentation
GTDB=""
if [ -z "!{EXTRACTED_DB}" ]
then 
  DATABASE=!{params.databases}/gtdb
  LOCK_FILE=${DATABASE}/lock.txt

  # Download plsdb database if necessary
  mkdir -p ${DATABASE}
  flock ${LOCK_FILE} concurrentDownload.sh --output=${DATABASE} \
	--link=!{DOWNLOAD_LINK} \
	--httpsCommand="wget -O gtdb.tar.gz !{DOWNLOAD_LINK} && tar xzvf gtdb.tar.gz && rm gtdb.tar.gz" \
	--s3FileCommand="s5cmd !{S5CMD_PARAMS} cp --concurrency !{task.cpus} !{DOWNLOAD_LINK} gtdb.tar.gz  && tar xzvf gtdb.tar.gz && rm gtdb.tar.gz " \
        --s3DirectoryCommand="s5cmd !{S5CMD_PARAMS} cp --concurrency !{task.cpus} !{DOWNLOAD_LINK} . " \
	--s5cmdAdditionalParams="!{S5CMD_PARAMS}" \
	--localCommand="tar -xzvf !{DOWNLOAD_LINK} " \
	--expectedMD5SUM=!{MD5SUM}

  GTDB=$(readlink -f ${DATABASE}/out/*)
else
  GTDB=!{EXTRACTED_DB}
fi

export GTDBTK_DATA_PATH=${GTDB}
gtdbtk classify_wf --batchfile input.tsv --out_dir output --cpus !{task.cpus} \
	--extension !{ending} !{GTDB_PARAMS}

# reformat gtdbtk output files
touch output/gtdbtk.bac120.summary.tsv
touch output/gtdbtk.ar122.summary.tsv
FILE_ID=$(mktemp -u XXXXXXXXXX)
FILE_BAC=chunk_${FILE_ID}_!{sample}_gtdbtk.bac120.summary.tsv
FILE_ARC=chunk_${FILE_ID}_!{sample}_gtdbtk.ar122.summary.tsv
FILE_COMB=!{sample}_gtdbtk_${FILE_ID}.tsv
FILE_UNCLASSIFIED=!{sample}_gtdbtk_unclassifed_${FILE_ID}.tsv

# Filter out unclassified
head -n 1 output/gtdbtk.bac120.summary.tsv > output/unclassified.tsv
grep "$(printf '\t')Unclassified$(printf '\t')" output/gtdbtk.bac120.summary.tsv >> output/unclassified.tsv || true
grep -v "$(printf '\t')Unclassified$(printf '\t')" output/gtdbtk.bac120.summary.tsv > output/classifiedBacteria.tsv || true

sed "s/^/SAMPLE\t/g" <(head -n 1 output/unclassified.tsv) > ${FILE_UNCLASSIFIED}
sed "s/^/!{sample}\t/g"  <(tail -n +2 output/unclassified.tsv) >> ${FILE_UNCLASSIFIED}

sed "s/^/SAMPLE\t/g" <(head -n 1 output/classifiedBacteria.tsv) > $FILE_BAC
sed "s/^/!{sample}\t/g"  <(tail -n +2 output/classifiedBacteria.tsv) >> $FILE_BAC

sed "s/^/SAMPLE\t/g" <(head -n 1 output/gtdbtk.ar122.summary.tsv) > $FILE_ARC
sed "s/^/!{sample}\t/g" <(tail -n +2 output/gtdbtk.ar122.summary.tsv) >> $FILE_ARC

GTDB_SUMMARY_TMP=gtdbtk_tmp.tsv
cat <(head -n 1 ${FILE_BAC}) <(head -n 1 ${FILE_ARC}) | sort | uniq | sed 's/^/DOMAIN\t/g' > $GTDB_SUMMARY_TMP
cat <(tail -n +2  ${FILE_ARC} | sed 's/^/ARCHAEA\t/g') <(tail -n +2  ${FILE_BAC} | sed 's/^/BACTERIA\t/g')  >> $GTDB_SUMMARY_TMP
paste -d$'\t' <(cut -f 3 $GTDB_SUMMARY_TMP | sed '1,1s/user_genome/BIN_ID/') $GTDB_SUMMARY_TMP > $FILE_COMB

# Copy trees to output
for tree in $(ls -1 output/classify/*.tree); do 
	NEW_TREE_NAME=$(basename $(echo ${tree} | sed "s/classify.tree/classify_${FILE_ID}.tree/g"));
        mv ${tree} ${NEW_TREE_NAME}
done
