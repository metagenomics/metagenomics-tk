readlink -f !{list_of_representatives} > list.txt 

additionalParams=" !{percentIdentity} !{covermParams} " 
covermCountParams=" --min-covered-fraction 0 "

coverm genome -t !{task.cpus} ${additionalParams} -b !{mapping} \
	--genome-fasta-list list.txt --methods mean --output-file !{sample}_mean.tsv && sed -i '1 s| Mean$||' !{sample}_mean.tsv
coverm genome -t !{task.cpus} ${additionalParams} -b !{mapping} \
	--genome-fasta-list list.txt --methods trimmed_mean --output-file !{sample}_trimmed_mean.tsv && sed -i '1 s| Trimmed Mean$||' !{sample}_trimmed_mean.tsv
coverm genome -t !{task.cpus} -b !{mapping} ${covermCountParams} \
	--genome-fasta-list list.txt --methods count  --output-file !{sample}_count.tsv  && sed -i '1 s| Read Count$||' !{sample}_count.tsv
coverm genome -t !{task.cpus} ${additionalParams} -b !{mapping} \
	--genome-fasta-list list.txt  --methods rpkm  --output-file !{sample}_rpkm.tsv  && sed -i '1 s| RPKM$||' !{sample}_rpkm.tsv
coverm genome -t !{task.cpus} ${additionalParams} -b !{mapping} \
	--genome-fasta-list list.txt --methods tpm --output-file !{sample}_tpm.tsv && sed -i '1 s| TPM$||' !{sample}_tpm.tsv
coverm genome -t !{task.cpus} ${additionalParams} -b !{mapping} \
	--genome-fasta-list list.txt --methods relative_abundance --output-file !{sample}_relative_abundance.tsv \
	&& sed -i '1 s| Relative Abundance (%)$||' !{sample}_relative_abundance.tsv
