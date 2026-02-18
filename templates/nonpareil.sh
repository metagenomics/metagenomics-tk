# Prepare fastq input for nonpareil by setting minimum sequence length and taking just one part of the paired end reads
FASTQ_MIN_LENGTH=24
zcat ${unpairedReads} | seqkit seq -j ${task.cpus} -m \${FASTQ_MIN_LENGTH} > unpaired.fq
zcat ${interleavedReads} | paste - - - - - - - -  \\
		      | cut -f 5-8 | tr "\\t" "\\n" | pigz --best --processes ${task.cpus} \\
		      | seqkit seq -j ${task.cpus} -m \${FASTQ_MIN_LENGTH} >> unpaired.fq

# Adjust the minimum number of required sequences if necessary
NUM_SEQS=\$(seqkit stats -T unpaired.fq | csvtk -t -T cut -f "num_seqs" | tail -n 1)
# Default is 10000
PERCENT_10="\$(( \${NUM_SEQS}*10/100 ))"
REQUIRED_SEQS="10000"
if [ \${PERCENT_10} -lt \${REQUIRED_SEQS} ]; then
  REQUIRED_SEQS="\${PERCENT_10}"
fi

# Run Nonpareil
nonpareil -t ${task.cpus} -R ${task.memory} -X \${REQUIRED_SEQS} -s unpaired.fq \\
	-T kmer -f fastq -b ${sample} ${params.steps.qc.nonpareil.additionalParams}

# Calculate nonpareil diversity index and other metrics
NPO=\$(readlink -f *.npo)
NONPAREIL_DIVERSITY_INDEX=${sample}_nonpareil_index.tsv
touch \${NONPAREIL_DIVERSITY_INDEX}
OUTPUT_PATH=\$(readlink -f \${NONPAREIL_DIVERSITY_INDEX})
R -e "library(Nonpareil);
      nps <- Nonpareil.set('\${NPO}'); 
      npsSummary <- summary(nps);
      npsSummaryIndex <- cbind(SAMPLE = rownames(npsSummary), npsSummary);
      write.table(npsSummaryIndex, file='\${OUTPUT_PATH}', quote=FALSE, sep='\\t', row.names=FALSE);"

mv *.pdf ${sample}_nonpareil_curves.pdf
