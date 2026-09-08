METASPADES_OUTPUT_DIR=output
ASSEMBLY_OUTPUT=\${METASPADES_OUTPUT_DIR}/contigs.fasta
ASSEMBLY_GRAPH_FASTG_OUTPUT=\${METASPADES_OUTPUT_DIR}/assembly_graph.fastg
ASSEMBLY_GRAPH_GFA_OUTPUT=\${METASPADES_OUTPUT_DIR}/assembly_graph_with_scaffolds.gfa
ASSEMBLY_GRAPH_PATHS_OUTPUT=\${METASPADES_OUTPUT_DIR}/contigs.paths

# run metaspades
spades.py -t ${task.cpus} --memory ${memory} \\
	--meta -o \${METASPADES_OUTPUT_DIR} --12 interleaved.fq.gz ${params.steps.assembly.metaspades.additionalParams}

ASSEMBLY_GZIPPED_OUTPUT=${sample}_contigs.fa.gz
HEADER_MAPPING_OUTPUT=${sample}_contigs_header_mapping.tsv

# The following function modifies the assembly fasta headers according to the pattern: SAMPLEID_SEQUENCECOUNTER_SEQUENCEHASH
transform.sh \${ASSEMBLY_OUTPUT} \${ASSEMBLY_GZIPPED_OUTPUT} \${HEADER_MAPPING_OUTPUT} ${sample} ${task.cpus}

# get basic contig stats 
paste -d\$'\\t' <(echo -e "SAMPLE\\n${sample}") <(seqkit stat -Ta \${ASSEMBLY_GZIPPED_OUTPUT}) > ${sample}_contigs_stats.tsv

cat > rename_paths.awk << 'AWKEOF'
NR == FNR {
    map[\$1] = \$2
    next
}
{
    name = \$0
    suffix = ""
    if (substr(name, length(name), 1) == "'") {
        suffix = "'"
        name = substr(name, 1, length(name) - 1)
    }
    if (name in map) {
        print map[name] suffix
    } else {
        print \$0
    }
}
AWKEOF

awk -f rename_paths.awk \${HEADER_MAPPING_OUTPUT} \${ASSEMBLY_GRAPH_PATHS_OUTPUT} > ${sample}_contigs_renamed.paths

# Export assembly graph
mv \${ASSEMBLY_GRAPH_GFA_OUTPUT} ${sample}_contigs.gfa
mv \${ASSEMBLY_GRAPH_PATHS_OUTPUT} ${sample}_contigs.paths

# transform assembly to assembly graph
maxKmer="default"
if [[ "${outputFastg}" == "TRUE" ]]; then
	# Maximum chosen Kmer
	maxKmer=\$(ls -1  \${METASPADES_OUTPUT_DIR}* | grep "^K" | sed 's/K//g' | sort -n | tail -n 1)
	mv \${ASSEMBLY_GRAPH_FASTG_OUTPUT} ${sample}_contigs.fastg
fi
