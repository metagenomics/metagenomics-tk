
# Get column with values and skip the column with names
for s in sample* ; do  
	cut -f 2- $s > ${s}_modified ; 
done

# Merge all files and add one column with names
paste -d$"\t" <(cut -f 1 sample1) *_modified >  abundance.tsv

