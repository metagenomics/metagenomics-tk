METASPADES_OUTPUT_DIR=output
ASSEMBLY_OUTPUT=${METASPADES_OUTPUT_DIR}/contigs.fasta
ASSEMBLY_GRAPH_OUTPUT=${METASPADES_OUTPUT_DIR}/assembly_graph.fastg

# run metaspades
spades.py -t !{task.cpus} --meta -o ${METASPADES_OUTPUT_DIR} --12 interleaved.fq.gz !{params.steps.assembly.metaspades.additionalParams}

ASSEMBLY_GZIPPED_OUTPUT=!{sample}_contigs.fa.gz

# The following function modifies the assembly fasta headers according to the pattern: SAMPLEID_SEQUENCECOUNTER_SEQUENCEHASH
transform.sh ${ASSEMBLY_OUTPUT} ${ASSEMBLY_GZIPPED_OUTPUT} !{sample} !{task.cpus}

# get basic contig stats 
paste -d$'\t' <(echo -e "SAMPLE\n!{sample}") <(seqkit stat -Ta ${ASSEMBLY_GZIPPED_OUTPUT}) > !{sample}_contigs_stats.tsv

# transform assembly to assembly graph
maxKmer="default"
if [[ "!{outputFastg}" == "TRUE" ]]; then
	# Maximum chosen Kmer
	maxKmer=$(ls -1  ${METASPADES_OUTPUT_DIR}* | grep "^K" | sed 's/K//g' | sort -n | tail -n 1)
	mv ${ASSEMBLY_GRAPH_OUTPUT} !{sample}_contigs.fastg
fi
