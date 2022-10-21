minimap2 !{params.steps.readMapping.minimap2.additionalParams.minimap2} -a !{representatives_fasta}  -t !{task.cpus} !{sample} - \
      | samtools view -@ !{task.cpus} -S -b - \
      | samtools sort -l 9 -@ !{task.cpus} - > !{sampleID}.bam


samtools index !{sampleID}.bam
