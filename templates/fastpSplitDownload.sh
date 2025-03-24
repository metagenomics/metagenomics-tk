# run fastp

set -o pipefail

s5cmd !{params.steps.qc.fastp.download.s5cmdParams} cat  --concurrency !{task.cpus} ${read1Url} 2> error1.log  > inputReads1.fq.gz
s5cmd !{params.steps.qc.fastp.download.s5cmdParams} cat  --concurrency !{task.cpus} ${read2Url} 2> error2.log  > inputReads2.fq.gz

fastp -i inputReads1.fq.gz \
      -I inputReads2.fq.gz \
      -o read1.fastp.fq.gz -O read2.fastp.fq.gz -w !{task.cpus} -h !{sample}_report.html \
         --unpaired1 !{sample}_unpaired.fastp.fq.gz --unpaired2 !{sample}_unpaired.qc.fq.gz !{params.steps.qc.fastp.additionalParams.fastp}

# This if statement solves issue https://github.com/pbelmann/meta-omics-toolkit/issues/166
if grep -q "reset by peer" error1.log error2.log; then
       echo "Network issue found. Exiting with exit code 1";
       exit 1 ;
else
       echo "No network issue found";
fi

# fix 'unexpected end of file' of unpaired reads gzip file
touch empty.txt
gzip empty.txt

UNPAIRED=!{sample}_unpaired.qc.fq.gz
cat empty.txt.gz >> ${UNPAIRED}

# create statistics for unpaired fastq files
paste -d$'\t' <(echo -e "SAMPLE\n!{sample}") <(seqkit stats -T ${UNPAIRED}) > !{sample}_unpaired_summary.tsv

# create interleaved fastq file for further analysis
paste <(zcat read1.fastp.fq.gz)  <(zcat read2.fastp.fq.gz) \
       | paste - - - - \
       | awk -v OFS="\n" -v FS="\t" '{print($1,$3,$5,$7,$2,$4,$6,$8)}' \
       | pigz --best --processes !{task.cpus} > !{sample}_interleaved.qc.fq.gz

# create tables of the fastp summary
cat fastp.json | jq -r  ' [.summary.before_filtering] | (map(keys) | add | unique) as $cols | map(. as $row | $cols | map($row[.])) as $rows | $cols, $rows[] | @tsv ' > fastp_summary_before_tmp.tsv
paste -d$'\t' <(echo -e "SAMPLE\n!{sample}") fastp_summary_before_tmp.tsv > !{sample}_fastp_summary_before.tsv

cat fastp.json | jq -r  ' [.summary.after_filtering] | (map(keys) | add | unique) as $cols | map(. as $row | $cols | map($row[.])) as $rows | $cols, $rows[] | @tsv ' > fastp_summary_after_tmp.tsv
paste -d$'\t' <(echo -e "SAMPLE\n!{sample}") fastp_summary_after_tmp.tsv > !{sample}_fastp_summary_after.tsv

mv fastp.json !{sample}_fastp.json
