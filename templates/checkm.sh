
# Prepare checkm patch, output directory and output file name
echo '{"dataRoot": "/.checkm", "remoteManifestURL": "https://data.ace.uq.edu.au/public/CheckM_databases/", "manifestType": "CheckM", "remoteManifestName": ".dmanifest", "localManifestName": ".dmanifest"}' > /tmp/DATA_CONFIG
mkdir out
FILE=$(mktemp !{sample}_checkm_XXXXXXXX.tsv)

# run suggested checkm commands
checkm tree !{params.steps.magAttributes.checkm.additionalParams.tree} --pplacer_threads !{task.cpus}  -t !{task.cpus} -x !{ending} . out &> tree.log
checkm tree_qa out &> tree_qa.log
checkm lineage_set !{params.steps.magAttributes.checkm.additionalParams.lineage_set} out out/marker &> lineage.log
checkm analyze -x !{ending} -t !{task.cpus} out/marker . out &> analyze.log
checkm qa !{params.steps.magAttributes.checkm.additionalParams.qa} --tab_table -t !{task.cpus} -f checkm.txt out/marker out  &> qa.log

# reformat output files according to magAttributes standard
echo -e "SAMPLE\tBIN_ID\tMarker lineage\t# genomes\t# markers\t# marker sets\t0\t1\t2\t3\t4\t5+\tCOMPLETENESS\tCONTAMINATION\tHETEROGENEITY" > $FILE
sed -i " 2,$ s/\t/.fa\t/"  checkm.txt
tail -n +2 checkm.txt | sed "s/^/!{sample}\t/g"  >> $FILE
