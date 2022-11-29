minimap2 !{params.steps.readMapping.minimap.additionalParams.minimap} -a !{index} -t !{task.cpus}  <(cat !{sample}) - \
      | samtools view -@ !{task.cpus} -S -b - \
      | samtools sort -l 9 -@ !{task.cpus} - > !{sampleID}.bam

samtools index !{sampleID}.bam
