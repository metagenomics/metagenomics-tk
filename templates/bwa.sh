gunzip -dc !{sample} | bwa mem -p -t !{task.cpus} !{representatives_fasta} - |  samtools sort -@ !{task.cpus} -o !{bin_shuffle_id}_!{ID}_bam.sorted; 
samtools index -@ !{task.cpus} !{bin_shuffle_id}_!{ID}_bam.sorted
