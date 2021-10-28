PLASMIDS_OUTPUT=!{sample}_plasmids.fasta.gz
scapp -g !{assemblyGraph} -k !{maxKmer} -p !{task.cpus} !{params.steps.plasmid.SCAPP.additionalParams} -b !{bam} -o .

# The following function modifies the assembly fasta headers according to the pattern: SAMPLEID_SEQUENCECOUNTER_SEQUENCEHASH
transform.sh "$(ls -1 *.confident_cycs.fasta)" ${PLASMIDS_OUTPUT} !{sample}

# get basic contig stats
csvtk concat --out-tabs -H <(csvtk transpose <(echo -e "SAMPLE\n!{sample}") ) <(csvtk transpose  <(seqkit stat -Ta ${PLASMIDS_OUTPUT}) ) \
	        | csvtk --tabs transpose > !{sample}_plasmids_stats.tsv
