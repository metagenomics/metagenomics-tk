# run megahit
megahit -t !{task.cpus} --12 interleaved.fq.gz  !{includeUnpairedReads} !{params.steps.assembly.megahit.additionalParams}

ASSEMBLY_OUTPUT=!{sample}_contigs.fa.gz
HEADER_MAPPING_OUTPUT=!{sample}_contigs_header_mapping.tsv
MEGAHIT_OUTPUT_DIR=megahit_out

# The following function modifies the assembly fasta headers according to the pattern: SAMPLEID_SEQUENCECOUNTER_SEQUENCEHASH
transform.sh ${MEGAHIT_OUTPUT_DIR}/final.contigs.fa ${ASSEMBLY_OUTPUT} ${HEADER_MAPPING_OUTPUT} !{sample} !{task.cpus}

# get basic contig stats 
paste -d$'\t' <(echo -e "SAMPLE\n!{sample}") <(seqkit stat -Ta ${ASSEMBLY_OUTPUT}) > !{sample}_contigs_stats.tsv

# transform assembly to assembly graph
maxKmer="default"
if [[ "!{convertToFastg}" == "TRUE" ]]; then
	# Maximum chosen Kmer
	maxKmer=$(cat ${MEGAHIT_OUTPUT_DIR}/options.json  | jq  '.k_max')

	# Use intermediate contigs of largest kmer for building the fastg file
	intermediateContigs="${MEGAHIT_OUTPUT_DIR}/intermediate_contigs/k${maxKmer}.contigs.fa"
	megahit_toolkit contig2fastg $maxKmer $intermediateContigs > !{sample}_contigs.fastg
fi

# Fix for ownership issue https://github.com/nextflow-io/nextflow/issues/4565
chmod -R  a+rw ${MEGAHIT_OUTPUT_DIR} 
