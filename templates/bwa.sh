gunzip -dc !{sample} | bwa mem !{MODE} -t !{task.cpus} !{representatives_fasta} - \
      | samtools view -@ !{task.cpus} -S -b - \
      | samtools sort -l 9 -@ !{task.cpus} - !{bin_shuffle_id}_!{ID} 
