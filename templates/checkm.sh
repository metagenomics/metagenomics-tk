
# Prepare checkm patch, output directory and output file name
echo '{"dataRoot": "/.checkm", "remoteManifestURL": "https://data.ace.uq.edu.au/public/CheckM_databases/", "manifestType": "CheckM", "remoteManifestName": ".dmanifest", "localManifestName": ".dmanifest"}' > /tmp/DATA_CONFIG
mkdir out
FILE=$(mktemp chunk_XXXXXXXXXX_!{sample}_checkm.txt)

# run suggested checkm commands
checkm tree --reduced_tree --pplacer_threads !{task.cpus}  -t !{task.cpus} -x !{ending} . out &> tree.log
checkm tree_qa out &> tree_qa.log
checkm lineage_set out out/marker &> lineage.log
checkm analyze -x !{ending} -t !{task.cpus} out/marker . out &> analyze.log
checkm qa --tab_table -t !{task.cpus} -f checkm.txt out/marker out  &> qa.log

# reformat output files
echo "SAMPLE\tBIN_ID\tMarker lineage\t# genomes\t# markers\t# marker sets\t0\t1\t2\t3\t4\t5+\tCOMPLETENESS\tCONTAMINATION\tHETEROGENEITY" > checkm_tmp.tsv
tail -n +2 checkm.txt | sed "s/^/!{sample}\t/g"  >> checkm_tmp.tsv
echo "PATH" > path.tsv
tail -n +2 checkm.txt | cut -f 1 | sed "s/$/!{ending}/g" | xargs -I {} readlink -f {} >> path.tsv
paste -d$'\t' path.tsv checkm_tmp.tsv > $FILE
