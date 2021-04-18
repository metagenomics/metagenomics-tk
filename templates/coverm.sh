OUT=!{ID}_!{sample}_out
mkdir $OUT
coverm genome -t !{task.cpus} --min-covered-fraction 0 -b !{mapping} --genome-fasta-list !{list_of_representatives} --methods mean --output-file $OUT/mean.tsv || true
coverm genome -t !{task.cpus} --min-covered-fraction 0 -b !{mapping} --genome-fasta-list !{list_of_representatives} --methods trimmed_mean --output-file $OUT/trimmed_mean.tsv || true
coverm genome -t !{task.cpus}  --min-covered-fraction 0 -b !{mapping} --genome-fasta-list !{list_of_representatives} --methods count  --output-file $OUT/count.tsv || true
coverm genome -t !{task.cpus} --min-covered-fraction 0 -b !{mapping} --genome-fasta-list !{list_of_representatives} --methods rpkm  --output-file $OUT/rpkm.tsv || true
coverm genome -t !{task.cpus} --min-covered-fraction 0 -b !{mapping} --genome-fasta-list !{list_of_representatives} --methods tpm --output-file $OUT/tpm.tsv || true
coverm genome -t !{task.cpus} --min-covered-fraction 10 -b !{mapping} --genome-fasta-list !{list_of_representatives} --methods mean --output-file $OUT/mean_mincov10.tsv || true
