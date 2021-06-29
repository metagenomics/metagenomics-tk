
# create input, output files and run default gtdbtk command
mkdir output
readlink -f !{bins} > bin.path
paste -d$'\t' bin.path <(for p in $(cat bin.path); do basename $p; done) > input.tsv

gtdbtk classify_wf --batchfile input.tsv --out_dir output --cpus !{task.cpus}  --extension !{ending}

# reformat gtdbtk output files
touch output/gtdbtk.bac120.summary.tsv
touch output/gtdbtk.ar122.summary.tsv
FILE_BAC=$(mktemp chunk_XXXXXXXXXX_!{sample}_gtdbtk.bac120.summary.tsv)
FILE_ARC=$(mktemp chunk_XXXXXXXXXX_!{sample}_gtdbtk.ar122.summary.tsv)
FILE_COMB=$(mktemp !{sample}_gtdbtk_XXXXXXX.tsv)

sed "s/^/SAMPLE\t/g" <(head -n 1 output/gtdbtk.bac120.summary.tsv) > $FILE_BAC
sed "s/^/!{sample}\t/g"  <(tail -n +2 output/gtdbtk.bac120.summary.tsv) >> $FILE_BAC

sed "s/^/SAMPLE\t/g" <(head -n 1 output/gtdbtk.ar122.summary.tsv) > $FILE_ARC
sed "s/^/!{sample}\t/g" <(tail -n +2 output/gtdbtk.ar122.summary.tsv) >> $FILE_ARC

cat <(head -n 1 ${FILE_BAC}) <(head -n 1 ${FILE_ARC}) | sort | uniq | sed 's/^/DOMAIN\t/g' > $FILE_COMB
cat <(tail -n +2  ${FILE_ARC} | sed 's/^/ARCHAEA\t/g') <(tail -n +2  ${FILE_BAC} | sed 's/^/BACTERIA\t/g')  >> $FILE_COMB
