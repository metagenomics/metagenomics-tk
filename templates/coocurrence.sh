
# Get BIN_ID, Domain and classification and column and replace semicolon with tab
cut -f 5 gtdb.input.tsv | sed 's/;/\t/g' \
	| sed '1 s/classification/DOMAIN\tPHYLUM\tCLASS\tORDER\tFAMILY\tGENUS\tSPECIES/g' > classification.tsv

cut -f 1 gtdb.input.tsv  | rev \
	| cut -f 2- -d '.' | rev > binId.tsv

paste -d$'\t' classification.tsv binId.tsv > gtdb.tsv

# Run the actual R coocurrence script
Rscript /cooccurrence.R -o $(pwd) -m abundance.tsv  -g gtdb.tsv -c !{task.cpus} !{params.steps.cooccurrence.additionalParams}
