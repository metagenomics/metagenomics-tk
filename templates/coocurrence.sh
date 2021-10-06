
# Get BIN_ID, Domain and classification and column and replace semicolon with tab
cut -f 1,2,5 gtdb.tsv \
	| sed 's/;/\t/g' \
	| sed '1 s/classification/DOMAIN\tPHYLUM\tCLASS\tORDER\tFAMILY\tGENUS\tSPECIES/g' > gtdb.fixed.tsv

# Run the actual R coocurrence script
Rscript /cooccurrence.R -o $(pwd) -m abundance.tsv  -g gtdb.fixed.tsv -c !{task.cpus} !{params.steps.cooccurrence.params}
