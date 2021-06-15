runMetaBat.sh !{contigs} !{bam}
TYPE=out
mkdir ${TYPE}
for bin in $(basename !{contigs})*/bin* ; do mv $bin ${TYPE}/!{sample}_$(basename ${bin}) ; done
for bin in ${TYPE}/*; do  BIN=$(basename $bin); grep ">" $bin | sed 's/>//g' | sed "s/^/${BIN}\t/g" >> bin_contig_mapping.tsv; done
head -n 1 *.depth.txt  | sed "s/$/\tBIN_ID\tSAMPLE/g" > !{sample}_bins_depth.tsv
join -t$'\t' -1 1 -2 2  <(tail -n +2  *.depth.txt | sort -k 1,1 ) <(sort -k 2,2 bin_contig_mapping.tsv) | sed "s/$/\t!{sample}/g" >> !{sample}_bins_depth.tsv
cd ${TYPE}
seqkit stat -Ta *.fa | sed '1s/^/SAMPLE\t/' | sed "1 ! s/^/!{sample}\t/" > ../!{sample}_bins_stats_tmp.tsv
cd ..

paste -d$'\t'  <(cut -f 2 !{sample}_bins_stats_tmp.tsv | sed 's/\.fa$//g' | tail -n +2 | sed '1s;^;BIN_ID\n;') !{sample}_bins_stats_tmp.tsv >  !{sample}_bins_stats_correct_ids.tsv

for id in  $(cut -f 6 !{sample}_bins_depth.tsv | tail -n +2 | sort | uniq); do  
	mean=$(awk -v bin=$id ' $6 == bin  { total += $3; count++ } END { print total/count }' !{sample}_bins_depth.tsv); echo -e "$id\t$mean" ;
done | sed '1s/^/BIN_ID\tCOVERAGE\n/' > avg_depth.tsv

join -t$'\t' --header -1 3 -2 1 !{sample}_bins_stats_correct_ids.tsv avg_depth.tsv > !{sample}_bins_stats.tsv
