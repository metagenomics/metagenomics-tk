PIGZ_COMPRESSION_THREADS=!{task.cpus}
zcat interleaved.fq.gz | paste  - - - - - - - - | tee >(cut -f 1-4 | tr "\t" "\n" | pigz --best --processes ${PIGZ_COMPRESSION_THREADS} > read1_1.fq.gz) | cut -f 5-8 | tr "\t" "\n" | pigz --best --processes ${PIGZ_COMPRESSION_THREADS} > read2_1.fq.gz
fastp -i read1_1.fq.gz -I read2_1.fq.gz -o read1_1.fastp.fq.gz -O read2_1.fastp.fq.gz -w !{task.cpus} -h !{sample}_report.html

cat fastp.json | jq -r  ' [.summary.before_filtering] | (map(keys) | add | unique) as $cols | map(. as $row | $cols | map($row[.])) as $rows | $cols, $rows[] | @tsv ' > fastp_summary.tsv

megahit -t !{task.cpus} -1  read1_1.fastp.fq.gz  -2 read2_1.fastp.fq.gz
paste <(zcat read1_1.fastp.fq.gz)  <(zcat read2_1.fastp.fq.gz) | paste - - - - | awk -v OFS="\n" -v FS="\t" '{print($1,$3,$5,$7,$2,$4,$6,$8)}'  | pigz --best --processes ${PIGZ_COMPRESSION_THREADS} > interleaved.fastp.fq.gz
TYPE="megahit"
mkdir !{sample}

mv megahit_out/final.contigs.fa !{sample}/

