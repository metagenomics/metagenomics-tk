    paste <(zcat read1.fq.gz)  <(zcat read2.fq.gz) \
	           | paste - - - - \
		          | awk -v OFS="\n" -v FS="\t" '{print($1,$3,$5,$7,$2,$4,$6,$8)}' \
			         | pigz --best --processes !{task.cpus} > !{sample}_interleaved.fq.gz

