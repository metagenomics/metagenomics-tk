bwa-mem2 mem !{params.steps.readMapping.bwa2.additionalParams.bwa2_mem} -p -t !{task.cpus} !{representatives} <(cat !{sample}) - \
      | samtools view -@ !{task.cpus} -S -b - \
      | samtools sort -l 9 -@ !{task.cpus} - > !{sampleID}.bam

samtools index !{sampleID}.bam
