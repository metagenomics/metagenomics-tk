mkdir refinement
select_representative.py -i !{genome_table} -d !{distance} -c !{cluster} -o refinement
cp refinement/representatives.tsv intermediate_clusters.tsv
