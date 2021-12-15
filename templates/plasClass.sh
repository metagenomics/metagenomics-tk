contigs="!{binID}_contigs.fa"
probabilities=!{sample}_!{binID}_probs.tsv
SEQUENCE_PROBABILITIES=sequence_length.tsv

pigz  -f -d -c !{assembly} > ${contigs}
classify_fasta.py -f ${contigs} -o ${SEQUENCE_PROBABILITIES} -p !{task.cpus} !{params.steps.plasmid.PlasClass.additionalParams}

# If sequences could be found then add additional information
if [ -s ${SEQUENCE_PROBABILITIES} ]; then
  echo -e "SEQUENCE\tPROBABILITY\tLENGTH\tSAMPLE\tBIN_ID" > ${probabilities}
  csvtk join --tabs --no-header-row  --fields "1" ${SEQUENCE_PROBABILITIES} \
	 <(seqkit fx2tab --length --only-id --name ${contigs} | sort -k 1,1) \
	| sed "s/$/\t!{sample}\t!{binID}/g" >> ${probabilities}
fi

# Filter the found sequences
FILTERED_SEQS=filtered_sequences.tsv
tail -n +2 ${probabilities} \
	| awk -v threshold=!{params.steps.plasmid.PlasClass.threshold} '($2+0 > threshold) {print $0}' - > ${FILTERED_SEQS} 

if [ -s ${FILTERED_SEQS} ]; then
	seqkit grep -f <(cut -f 1 ${FILTERED_SEQS}) !{binID}_contigs.fa -o !{sample}_!{binID}_plasmids.fasta.gz
fi
