
# Get BIN_ID, Domain and classification and column and replace semicolon with tab
cut -f 5 gtdb.input.tsv | sed 's/;/\t/g' \
	| sed '1 s/classification/DOMAIN\tPHYLUM\tCLASS\tORDER\tFAMILY\tGENUS\tSPECIES/g' > classification.tsv

cut -f 1 gtdb.input.tsv  | rev \
	| cut -f 2- -d '.' | rev > binId.tsv

paste -d$'\t' classification.tsv binId.tsv > gtdb.tsv

# Run the actual R coocurrence script
Rscript /cooccurrence.R create -o $(pwd) -m abundance.tsv  -g gtdb.tsv -c !{task.cpus} !{params.steps.cooccurrence.inference.additionalParams}

# Prepend row number to every line which should allow to batch the computation
awk -v batches=!{batchSize} \
	'{ printf("%0.0f\t",(NR/(batches+.001))+.5) ;print $1"\t"$2 }' <(tail -n +2 edges.tsv | sort -t$'\t' -k 1,1) \
	| sed '1s;^;IDX\tV1\tV2\n;' > edges_index.tsv
