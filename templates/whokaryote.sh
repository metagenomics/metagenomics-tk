gunzip -c !{gff} > unzipped.gff

whokaryote.py --contigs !{fasta} --threads !{task.cpus} --gff unzipped.gff --outdir output !{params.steps.annotation.whokaryote.additionalParams} 

mv output/featuretable.csv !{binID}_featuretable.csv 
mv output/featuretable_predictions_T.tsv !{binID}_featuretable_predictions_T.tsv
mv output/eukaryote_contig_headers.txt !{binID}_eukaryote_contig_headers.txt
mv output/prokaryote_contig_headers.txt !{binID}_prokaryote_contig_headers.txt
