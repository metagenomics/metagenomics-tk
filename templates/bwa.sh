bwa mem !{params.steps.readMapping.bwa.additionalParams.bwa_mem} -p -t !{task.cpus} !{representatives} <(cat !{sample}) - \
      | samtools view -@ !{task.cpus} -S -b - \
      | samtools sort -l 9 -@ !{task.cpus} - > !{sampleID}.bam

samtools index !{sampleID}.bam
