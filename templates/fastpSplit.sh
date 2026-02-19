# run fastp
fastp -i read1.fq.gz -I read2.fq.gz ${reportOnly} \\
	-w ${task.cpus} -h ${sample}_report.html \\
	--unpaired1 ${sample}_tmp_unpaired.qc.fq.gz --unpaired2 ${sample}_tmp_unpaired.qc.fq.gz ${params.steps.qc.fastp.additionalParams.fastp}

# fix 'unexpected end of file' of unpaired reads gzip file
touch empty.txt
gzip empty.txt

UNPAIRED=${sample}_tmp_unpaired.qc.fq.gz
cat empty.txt.gz >> \${UNPAIRED}

# create statistics for unpaired fastq files
paste -d\$'\\t' <(echo -e "SAMPLE\\n${sample}") <(seqkit stats -T \${UNPAIRED}) > ${sample}_unpaired_summary.tsv

# create interleaved fastq file for further analysis
paste <(zcat read1.fastp.fq.gz)  <(zcat read2.fastp.fq.gz) \\
       | paste - - - - \\
       | awk -v OFS="\\n" -v FS="\\t" '{print(\$1,\$3,\$5,\$7,\$2,\$4,\$6,\$8)}' \\
       | pigz --best --processes ${task.cpus} > ${sample}_tmp_interleaved.qc.fq.gz

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
