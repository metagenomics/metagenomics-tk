METASPADES_OUTPUT_DIR=output
ASSEMBLY_OUTPUT=${METASPADES_OUTPUT_DIR}/contigs.fasta
ASSEMBLY_GRAPH_OUTPUT=${METASPADES_OUTPUT_DIR}/assembly_graph.fastg

# run metaspades
MEMORY=$(echo !{task.memory} | awk '{val=$1; unit=$2; if (unit == "TB") printf "%.0f\n", val*1000; else if (unit == "GB") printf "%.0f\n", val}')
spades.py -t !{task.cpus} --memory ${MEMORY} --meta -o ${METASPADES_OUTPUT_DIR} --12 interleaved.fq.gz --nanopore ontReads.fq.gz !{params.steps.assemblyHybrid.metaspades.additionalParams}

ASSEMBLY_GZIPPED_OUTPUT=!{sample}_contigs.fa.gz
HEADER_MAPPING_OUTPUT=!{sample}_contigs_header_mapping.tsv

# The following function modifies the assembly fasta headers according to the pattern: SAMPLEID_SEQUENCECOUNTER_SEQUENCEHASH
transform.sh ${ASSEMBLY_OUTPUT} ${ASSEMBLY_GZIPPED_OUTPUT} ${HEADER_MAPPING_OUTPUT} !{sample} !{task.cpus}

# get basic contig stats 
paste -d$'\t' <(echo -e "SAMPLE\n!{sample}") <(seqkit stat -Ta ${ASSEMBLY_GZIPPED_OUTPUT}) > !{sample}_contigs_stats.tsv

# transform assembly to assembly graph
maxKmer="default"
if [[ "!{outputFastg}" == "TRUE" ]]; then
	# Maximum chosen Kmer
	maxKmer=$(ls -1  ${METASPADES_OUTPUT_DIR}* | grep "^K" | sed 's/K//g' | sort -n | tail -n 1)
	mv ${ASSEMBLY_GRAPH_OUTPUT} !{sample}_contigs.fastg
fi
