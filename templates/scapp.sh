PLASMIDS_OUTPUT=!{sample}_plasmids.fasta.gz

# Make sure that exit code 1 is also accepted due to known SCAPP bug:
# https://github.com/Shamir-Lab/SCAPP/issues/26
trap 'if [[ $? == 1 ]]; then echo OK; exit 0; fi' EXIT

scapp -g !{assemblyGraph} -k !{maxKmer} -p !{task.cpus} !{params.steps.plasmid.SCAPP.additionalParams.SCAPP} -b !{bam} -o .

if [ -s *.confident_cycs.fasta ]; then
  # The following function modifies the assembly fasta headers according to the pattern: SAMPLEID_SEQUENCECOUNTER_SEQUENCEHASH
  transform.sh "$(ls -1 *.confident_cycs.fasta)" ${PLASMIDS_OUTPUT} !{sample} !{task.cpus}

  PLASMID_STATS=!{sample}_plasmids_stats.tsv
  # get basic contig stats
  csvtk concat --out-tabs -H <(csvtk transpose <(echo -e "SAMPLE\n!{sample}") ) <(csvtk transpose  <(seqkit stat -Ta ${PLASMIDS_OUTPUT}) ) \
	        | csvtk --tabs transpose | sed 's/"//g' > ${PLASMID_STATS}

  echo -e "SAMPLE\tID\tLENGTH\tGC" > ${PLASMID_STATS}
  seqkit fx2tab -l -g -n -i  ${PLASMIDS_OUTPUT} \
	| sed "s/^/!{sample}\t/g"  >> ${PLASMID_STATS} 
fi
