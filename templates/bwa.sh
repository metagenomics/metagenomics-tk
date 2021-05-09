gunzip -dc !{sample} | bwa mem -p -t 20 !{representatives_fasta} - |  samtools sort -@ 20 -o !{bin_shuffle_id}_!{ID}_bam.sorted; 
CONTIGS=$(grep ">" !{representatives_fasta} | tr -d '>' | tr '\n' ' ')
samtools index -@ 20 !{bin_shuffle_id}_!{ID}_bam.sorted
