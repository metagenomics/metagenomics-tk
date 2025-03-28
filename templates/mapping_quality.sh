base=!{sample}_flagstat
out=${base}.tsv
echo -e "!{sample}\t!{sample}\tSAMPLE" > $out
samtools flagstat -@ !{task.cpus} -O tsv !{bam} >> $out

cut -f 3 $out | paste -s -d '\t'  > ${base}_failed.tsv
cut -f 2 $out | paste -s -d '\t' >> ${base}_failed.tsv
cut -f 3 $out | paste -s -d '\t' > ${base}_passed.tsv
cut -f 1 $out | paste -s -d '\t' >> ${base}_passed.tsv

