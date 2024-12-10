annotatedcontigs2json -bins-dir bins \
	-sample-bam-files  !{mapping} \
	-fasta !{contigs} \
	-json-gz !{sample}.contigs.json.gz \
	-sample-names !{sample}
