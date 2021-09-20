
# Get BIN_ID, Domain and classification and column and replace semicolon with tab
cut -f 1,2,5 gtdb.tsv \ 
	| sed 's/;/\t/g' \ 
	| sed '1 s/classification/DOMAIN\tPHYLUM\tCLASS\tORDER\tFAMILY\tGENUS\tSPECIES/g' > gtdb.fixed.tsv

# Run the actual R coocurrence script
Rscript /cooccurrence.R -o $(pwd) -m abundance.tsv  -g gtdb.fixed.tsv -c !{task.cpus} -t 0.1 -z 80 -o "correlation"
