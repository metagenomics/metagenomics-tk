PIGZ_COMPRESSION_THREADS=!{task.cpus}
fastp -i read1.fq.gz -I read2.fq.gz -o read1.fastp.fq.gz -O read2.fastp.fq.gz -w !{task.cpus} -h !{sample}_report.html
cat fastp.json | jq -r  ' [.summary.before_filtering] | (map(keys) | add | unique) as $cols | map(. as $row | $cols | map($row[.])) as $rows | $cols, $rows[] | @tsv ' > fastp_summary.tsv

megahit -t !{task.cpus} -1  read1.fastp.fq.gz  -2 read2.fastp.fq.gz
TYPE="megahit"
paste <(zcat read1.fastp.fq.gz)  <(zcat read2.fastp.fq.gz) | paste - - - - | awk -v OFS="\n" -v FS="\t" '{print($1,$3,$5,$7,$2,$4,$6,$8)}'  | pigz --best --processes ${PIGZ_COMPRESSION_THREADS} > interleaved.fastp.fq.gz
CONTIG_NAME=!{sample}_$(echo $RANDOM)_$(echo $RANDOM)
cat megahit_out/final.contigs.fa | seqkit replace  -p '(^\w*)' -r "\${1}_$CONTIG_NAME" > final.contigs.fa
