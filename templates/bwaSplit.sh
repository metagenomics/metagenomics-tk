bwa mem -t !{task.cpus} !{representatives_fasta} <(gunzip -dc !{read1}) <(gunzip -dc !{read2}) | \
	samtools view -@ !{task.cpus} -S -b - | samtools sort -l 9 -@ !{task.cpus} - !{bin_shuffle_id} 

samtools index !{bin_shuffle_id}.bam
