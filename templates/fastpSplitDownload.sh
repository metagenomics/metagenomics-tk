# run fastp

set -o pipefail

if [[ "\${isInterleaved}" == "true" ]]; then
    s5cmd ${params.steps.qc.fastp.download.s5cmdParams} cat  --concurrency ${task.cpus} \${read1Url} 2> error1.log  \\
       | zcat | paste - - - - - - - -  | tee >(cut -f 1-4 | tr "\t" "\n" | pigz --best --processes ${task.cpus} > inputReads1.fq.gz)  \\
       | cut -f 5-8 | tr "\t" "\n" | pigz --best --processes ${task.cpus} > inputReads2.fq.gz
else
    s5cmd ${params.steps.qc.fastp.download.s5cmdParams} cat  --concurrency ${task.cpus} \${read1Url} 2> error1.log  > inputReads1.fq.gz
    s5cmd ${params.steps.qc.fastp.download.s5cmdParams} cat  --concurrency ${task.cpus} \${read2Url} 2> error2.log  > inputReads2.fq.gz
fi

fastp -i inputReads1.fq.gz \\
      -I inputReads2.fq.gz \\
      --stdout \\
      -w ${task.cpus} \\
      -h ${sample}_report.html \\
       --unpaired1 ${sample}_tmp_unpaired.qc.fq.gz --unpaired2 ${sample}_tmp_unpaired.qc.fq.gz ${params.steps.qc.fastp.additionalParams.fastp} \\
	| pigz --best --processes ${task.cpus} > ${sample}_tmp_interleaved.qc.fq.gz

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

UNPAIRED=${sample}_tmp_unpaired.qc.fq.gz
cat empty.txt.gz >> \${UNPAIRED}

# create statistics for unpaired fastq files
paste -d\$'\\t' <(echo -e "SAMPLE\\n${sample}") <(seqkit stats -T \${UNPAIRED}) > ${sample}_unpaired_summary.tsv

# create tables of the fastp summary
cat fastp.json | jq -r  ' [.summary.before_filtering] | (map(keys) | add | unique) as \$cols | map(. as \$row | \$cols | map(\$row[.])) as \$rows | \$cols, \$rows[] | @tsv ' > fastp_summary_before_tmp.tsv
paste -d\$'\\t' <(echo -e "SAMPLE\\n${sample}") fastp_summary_before_tmp.tsv > ${sample}_fastp_summary_before.tsv

cat fastp.json | jq -r  ' [.summary.after_filtering] | (map(keys) | add | unique) as \$cols | map(. as \$row | \$cols | map(\$row[.])) as \$rows | \$cols, \$rows[] | @tsv ' > fastp_summary_after_tmp.tsv
paste -d\$'\\t' <(echo -e "SAMPLE\\n${sample}") fastp_summary_after_tmp.tsv > ${sample}_fastp_summary_after.tsv

mv fastp.json ${sample}_fastp.json

# Make sure that resulting fastq files are not empty 
TOTAL_READS_AFTER=\$(csvtk -T -t cut -f total_reads ${sample}_fastp_summary_after.tsv | tail -n 1)
if [ "\$TOTAL_READS_AFTER" -ne 0 ]; then
  mv ${sample}_tmp_interleaved.qc.fq.gz ${sample}_interleaved.qc.fq.gz
  mv ${sample}_tmp_unpaired.qc.fq.gz ${sample}_unpaired.qc.fq.gz
fi
