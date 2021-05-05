fastp -i ${genomeReads1} -I ${genomeReads2} -o $fq_1_paired -O $fq_2_paired -w ${task.cpus} -h ${sample}_fastp_report.html -j ${sample}_fastp_report.json
cat fastp.json | jq -r  ' [.summary.before_filtering] | (map(keys) | add | unique) as $cols | map(. as $row | $cols | map($row[.])) as $rows | $cols, $rows[] | @tsv ' > fastp_summary.tsv
