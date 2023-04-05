contigs="!{binID}_contigs.fa"
PLASCLASS_OUT=out.tsv
SEQUENCE_PROBABILITIES=sequence_length.tsv
FINAL_OUTPUT=!{binID}_plasclass.tsv
TMP_OUTPUT=!{binID}_plasclass_tmp.tsv

pigz  -f -d -c !{assembly} > ${contigs}
classify_fasta.py -f ${contigs} -o ${SEQUENCE_PROBABILITIES} -p !{task.cpus} !{params.steps.plasmid.PlasClass.additionalParams}

# If sequences could be found then add additional information
if [ -s ${SEQUENCE_PROBABILITIES} ]; then
  HEADER="CONTIG\tPROBABILITY\tLENGTH\tSAMPLE\tBIN_ID"
  echo -e ${HEADER} > ${PLASCLASS_OUT}
  csvtk join --tabs --no-header-row  --fields "1" ${SEQUENCE_PROBABILITIES} \
	 <(seqkit fx2tab --length --only-id --name ${contigs} | sort -k 1,1) \
	| sed "s/$/\t!{sample}\t!{binID}/g" >> ${PLASCLASS_OUT}

  echo -e ${HEADER} > ${TMP_OUTPUT}

  # Filter the found sequences
  tail -n +2 ${PLASCLASS_OUT} \
        | awk -v threshold=!{params.steps.plasmid.PlasClass.threshold} '($2+0 > threshold) {print $0}' - >> ${TMP_OUTPUT}
  
  if [ $(wc -l < ${TMP_OUTPUT}) -gt 1 ]; then
         mv ${TMP_OUTPUT} ${FINAL_OUTPUT}
  fi
fi

