OUT=!{sample}_out
mkdir $OUT

readlink -f !{list_of_representatives} > list.txt 

coverm genome -t !{task.cpus} !{params.steps.readMapping.bwa.additionalParams.coverm} -b !{mapping} \
	--genome-fasta-list list.txt --methods mean --output-file $OUT/mean.tsv && sed -i '1 s| Mean$||' $OUT/mean.tsv  || true
coverm genome -t !{task.cpus} !{params.steps.readMapping.bwa.additionalParams.coverm} -b !{mapping} \
	--genome-fasta-list list.txt --methods trimmed_mean --output-file $OUT/trimmed_mean.tsv && sed -i '1 s| Trimmed Mean$||' $OUT/trimmed_mean.tsv || true
coverm genome -t !{task.cpus} !{params.steps.readMapping.bwa.additionalParams.coverm} -b !{mapping} \
	--genome-fasta-list list.txt --methods count  --output-file $OUT/count.tsv  && sed -i '1 s| Read Count$||' $OUT/count.tsv || true
coverm genome -t !{task.cpus} !{params.steps.readMapping.bwa.additionalParams.coverm} -b !{mapping} \
	--genome-fasta-list  list.txt  --methods rpkm  --output-file $OUT/rpkm.tsv  && sed -i '1 s| RPKM$||' $OUT/rpkm.tsv || true
coverm genome -t !{task.cpus} !{params.steps.readMapping.bwa.additionalParams.coverm} -b !{mapping} \
	--genome-fasta-list list.txt --methods tpm --output-file $OUT/tpm.tsv && sed -i '1 s| TPM$||' $OUT/tpm.tsv || true

