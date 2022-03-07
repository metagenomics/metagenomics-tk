gunzip -dc !{sample} | bwa mem !{MODE} !{params.steps.readMapping.bwa.additionalParams.bwa_mem} -t !{task.cpus} !{representatives_fasta} - \
      | samtools view -@ !{task.cpus} -S -b - \
      | samtools sort -l 9 -@ !{task.cpus} - > !{sampleID}.bam

samtools index !{sampleID}.bam
