gunzip -dc !{sample} | bwa mem !{params.steps.readMapping.bwa.additionalParams.bwa_mem} -p -t !{task.cpus} !{representatives_fasta} - \
	| samtools view -@ !{task.cpus} -S -b - \
	| samtools sort -l 9 -@ !{task.cpus} - !{bin_shuffle_id} 

samtools index !{bin_shuffle_id}.bam
