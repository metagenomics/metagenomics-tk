#!/bin/bash
# Create a script that creates for every bin a seperate file 

mkdir tmp

cat *.fas > tmp/all.fas

# Prepare fasta header sequence
seqkit fx2tab tmp/all.fas | cut -d$'\t' -f 1 > tmp/header.tsv
sed -i 's/|||||/\t/g' tmp/header.tsv

# Prepare table for  
paste -d$'\t' <(seqkit fx2tab tmp/all.fas | cut -d$'\t' -f 1) tmp/header.tsv > tmp/header_bin_id.tsv

sed -i '1s/^/CONTIG_CURRENT\tCONTIG_NEW\tCONTIG_OLD\tBIN_ID\n/' tmp/header_bin_id.tsv
csvtk -T -t cut -f CONTIG_CURRENT,CONTIG_NEW tmp/header_bin_id.tsv  > tmp/header_replacement.tsv
seqkit replace -p '^(\S+)' -r '{kv}' -k tmp/header_replacement.tsv tmp/all.fas  > tmp/all_updated_header.fas

csvtk cut -t -T -f BIN_ID tmp/header_bin_id.tsv | tail -n +2 \
| sort | uniq | xargs -P !{task.cpus} -I {} splitFas.sh "{}" "!{sample}_!{type}.!{dbType}.fas"
