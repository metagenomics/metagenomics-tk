
# create input, output files and run default gtdbtk command
mkdir output
ls -1 !{bins} > binNames.tsv
ls -1 !{bins} | xargs -I {} readlink -f {} > bin.path
paste -d$'\t' bin.path <(for p in $(cat bin.path); do basename $p; done) > input.tsv

gtdb_download.sh "!{EXTRACTED_DB}" "!{DOWNLOAD_LINK}" "!{S5CMD_PARAMS}" "!{task.cpus}" "!{params.polished.databases}" "!{MD5SUM}" "!{S3_gtdb_ACCESS}" "!{S3_gtdb_SECRET}" || exit 1
GTDB=$(cat gtdbPath.txt)

export GTDBTK_DATA_PATH=${GTDB}
gtdbtk classify_wf --batchfile input.tsv --out_dir output --cpus !{task.cpus} \
	--mash_db genomes.msh --extension !{ending} !{GTDB_PARAMS}

# reformat gtdbtk output files
touch output/gtdbtk.bac120.summary.tsv
touch output/gtdbtk.ar53.summary.tsv
FILE_ID=!{FILE_TYPE}!{chunkId}
FILE_BAC=chunk_${FILE_ID}_!{sample}_gtdbtk.bac120.summary.tsv
FILE_BAC_TMP=chunk_${FILE_ID}_!{sample}_gtdbtk.bac120.summary.tmp.tsv
FILE_ARC=chunk_${FILE_ID}_!{sample}_gtdbtk.ar53.summary.tsv
FILE_ARC_TMP=chunk_${FILE_ID}_!{sample}_gtdbtk.ar53.summary.tmp.tsv
FILE_COMB_TMP=chunk_${FILE_ID}_!{sample}_gtdbtk_combined.tmp.tsv
FILE_COMB=chunk_${FILE_ID}_!{sample}_gtdbtk_combined.tsv
FILE_SUMMARY_COMBINED=chunk_${FILE_ID}_!{sample}_gtdbtk_summary_combined.tsv
FILE_SUMMARY_RAW_COMBINED=chunk_${FILE_ID}_!{sample}_gtdbtk_summary_raw_combined.tsv
FILE_SUMMARY_RAW_COMBINED_TMP=chunk_${FILE_ID}_!{sample}_gtdbtk_summary_raw_combined_tmp.tsv
FILE_UNCLASSIFIED=chunk_${FILE_ID}_!{sample}_gtdbtk_unclassified.tsv
FILE_UNCLASSIFIED_TMP=chunk_${FILE_ID}_!{sample}_gtdbtk_unclassified.tmp.tsv
MISSING_OUTPUT=chunk_${FILE_ID}_!{sample}_missing_bins.tsv
MISSING_OUTPUT_TMP=chunk_${FILE_ID}_!{sample}_missing_bins.tmp.tsv

# Filter out unclassified
head -n 1 output/gtdbtk.bac120.summary.tsv > output/unclassified.tsv
grep "$(printf '\t')Unclassified$(printf '\t')" output/gtdbtk.bac120.summary.tsv >> output/unclassified.tsv || true
grep "$(printf '\t')Unclassified Bacteria$(printf '\t')" output/gtdbtk.bac120.summary.tsv >> output/unclassified.tsv || true
grep "$(printf '\t')Unclassified Archaea$(printf '\t')" output/gtdbtk.ar53.summary.tsv >> output/unclassified.tsv || true

grep -v "$(printf '\t')Unclassified$(printf '\t')" output/gtdbtk.bac120.summary.tsv \
       | grep -v "$(printf '\t')Unclassified Bacteria$(printf '\t')" > output/classifiedBacteria.tsv || true

grep -v "$(printf '\t')Unclassified Archaea$(printf '\t')" output/gtdbtk.ar53.summary.tsv > output/classifiedArchaea.tsv || true

sed "s/^/SAMPLE\t/g" <(head -n 1 output/unclassified.tsv) > ${FILE_UNCLASSIFIED_TMP}
sed "s/^/!{sample}\t/g"  <(tail -n +2 output/unclassified.tsv) >> ${FILE_UNCLASSIFIED_TMP}

sed "s/^/SAMPLE\t/g" <(head -n 1 output/classifiedBacteria.tsv) > $FILE_BAC_TMP
sed "s/^/!{sample}\t/g"  <(tail -n +2 output/classifiedBacteria.tsv) >> $FILE_BAC_TMP

sed "s/^/SAMPLE\t/g" <(head -n 1 output/classifiedArchaea.tsv) > $FILE_ARC_TMP
sed "s/^/!{sample}\t/g" <(tail -n +2 output/classifiedArchaea.tsv) >> $FILE_ARC_TMP

GTDB_SUMMARY_TMP=gtdbtk_tmp.tsv
cat <(head -n 1 ${FILE_BAC_TMP}) <(head -n 1 ${FILE_ARC_TMP}) | sort | uniq | sed 's/^/DOMAIN\t/g' > $GTDB_SUMMARY_TMP
cat <(tail -n +2  ${FILE_ARC_TMP} | sed 's/^/ARCHAEA\t/g') <(tail -n +2  ${FILE_BAC_TMP} | sed 's/^/BACTERIA\t/g')  >> $GTDB_SUMMARY_TMP
paste -d$'\t' <(cut -f 3 $GTDB_SUMMARY_TMP | sed '1,1s/user_genome/BIN_ID/') $GTDB_SUMMARY_TMP > $FILE_COMB_TMP

# Since EMGB importer doesn't work with additional columns, we also export the raw version of the summary (with the sample column).
# Once EMGB support additional columns and only reads data by reading column names, the following block is no longer necessary.
GTDB_SUMMARY_RAW_COMBINED=gtdbtk_raw_tmp.tsv
cat <(head -n 1 ${FILE_BAC_TMP}) <(head -n 1 ${FILE_ARC_TMP}) | sort | uniq > ${FILE_SUMMARY_RAW_COMBINED_TMP}
cat <(tail -n +2  ${FILE_ARC_TMP}) <(tail -n +2  ${FILE_BAC_TMP})  >> $FILE_SUMMARY_RAW_COMBINED_TMP
if [[ $(wc -l <${FILE_SUMMARY_RAW_COMBINED_TMP}) -ge 2 ]]; then
	mv $FILE_SUMMARY_RAW_COMBINED_TMP $FILE_SUMMARY_RAW_COMBINED
fi

if [[ $(wc -l <${GTDB_SUMMARY_TMP}) -ge 2 ]]; then
	        mv ${GTDB_SUMMARY_TMP} ${FILE_SUMMARY_COMBINED}
fi

if [[ $(wc -l <${FILE_COMB_TMP}) -ge 2 ]]; then
	        mv ${FILE_COMB_TMP} ${FILE_COMB}
fi

if [[ $(wc -l <${FILE_BAC_TMP}) -ge 2 ]]; then
	        mv ${FILE_BAC_TMP} ${FILE_BAC}
fi

if [[ $(wc -l <${FILE_ARC_TMP}) -ge 2 ]]; then
	        mv ${FILE_ARC_TMP} ${FILE_ARC}
fi

if [[ $(wc -l <${FILE_UNCLASSIFIED_TMP}) -ge 2 ]]; then
	        mv ${FILE_UNCLASSIFIED_TMP} ${FILE_UNCLASSIFIED}
fi

# Copy trees to output
for tree in $(find output/classify/ -name "*.tree"); do 
	NEW_TREE_NAME=$(basename $(echo ${tree} | sed "s/classify.tree/classify_${FILE_ID}.tree/g"));
        mv ${tree} ${NEW_TREE_NAME}
done

FOUND_TMP=found.tsv
if [[ -f "${FILE_COMB}" ]]; then
	 cp "${FILE_COMB}"  ${FOUND_TMP}
else
	 touch ${FOUND_TMP}
fi

echo -e "SAMPLE\tBIN_ID" > ${MISSING_OUTPUT_TMP}
cat binNames.tsv <(cut -f 1 ${FOUND_TMP} | tail -n +2) \
	| sort | uniq -c | tr -s ' ' | grep " 1 "  | cut -f 3 -d ' ' | sed "s/^/!{sample}\t/g"  >> ${MISSING_OUTPUT_TMP}

if [[ $(wc -l <${MISSING_OUTPUT_TMP}) -ge 2 ]]; then
	        mv  ${MISSING_OUTPUT_TMP} ${MISSING_OUTPUT}
fi

# Fix for ownership issue https://github.com/nextflow-io/nextflow/issues/4565
chmod a+rw -R output
