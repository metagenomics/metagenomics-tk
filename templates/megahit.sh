# run megahit
megahit -t !{task.cpus} --12  reads.fq.gz

# make contig identifiers unique
CONTIG_NAME=!{sample}_$(echo $RANDOM)_$(echo $RANDOM)
cat megahit_out/final.contigs.fa \
       	| seqkit replace  -p '(^\w*)' -r "\${1}_$CONTIG_NAME" \
       	| pigz --best --processes !{task.cpus}  > !{sample}_contigs.fa.gz

# get basic contig stats 
paste -d$'\t' <(echo -e "SAMPLE\n!{sample}") <(seqkit stat -Ta  !{sample}_final.contigs.fa.gz) > !{sample}_contigs_stats.tsv
