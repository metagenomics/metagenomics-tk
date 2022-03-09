# run fastp
fastp   --stdin -i <(s5cmd !{params.steps.qc.fastp.download.s5cmdParams} cat ${read1Url} | zcat) \
      -I <(s5cmd !{params.steps.qc.fastp.download.s5cmdParams} cat ${read2Url} | zcat) \
      -o read1.fastp.fq.gz -O read2.fastp.fq.gz -w !{task.cpus} -h !{sample}_report.html \
      --unpaired1 !{sample}_unpaired.fastp.fq.gz --unpaired2 !{sample}_unpaired.qc.fq.gz !{params.steps.qc.fastp.additionalParams}

# fix 'unexpected end of file' of unpaired reads gzip file
touch empty.txt
gzip empty.txt
cat empty.txt.gz >> !{sample}_unpaired.qc.fq.gz

# create interleaved fastq file for further analysis
paste <(zcat read1.fastp.fq.gz)  <(zcat read2.fastp.fq.gz) \
       | paste - - - - \
       | awk -v OFS="\n" -v FS="\t" '{print($1,$3,$5,$7,$2,$4,$6,$8)}' \
       | pigz --best --processes !{task.cpus} > !{sample}_interleaved.qc.fq.gz

# create tables of the fastp summary
cat fastp.json | jq -r  ' [.summary.before_filtering] | (map(keys) | add | unique) as $cols | map(. as $row | $cols | map($row[.])) as $rows | $cols, $rows[] | @tsv ' > fastp_summary_before_tmp.tsv
paste -d$'\t' <(echo -e "SAMPLE\n!{sample}") fastp_summary_before_tmp.tsv > fastp_summary_before.tsv

cat fastp.json | jq -r  ' [.summary.after_filtering] | (map(keys) | add | unique) as $cols | map(. as $row | $cols | map($row[.])) as $rows | $cols, $rows[] | @tsv ' > fastp_summary_after_tmp.tsv
paste -d$'\t' <(echo -e "SAMPLE\n!{sample}") fastp_summary_after_tmp.tsv > fastp_summary_after.tsv
