METACAT_INDEX_PARAMS="${params.steps.binning.metacat.additionalParams.index}"
METACAT_SEED_PARAMS="${params.steps.binning.metacat.additionalParams.seed}"
METACAT_COVERAGE_PARAMS="${params.steps.binning.metacat.additionalParams.coverage}"
METACAT_CLUSTER_PARAMS="${params.steps.binning.metacat.additionalParams.cluster}"

TEMP_DIR=\$(basename \$(mktemp))
mkdir \${TEMP_DIR}

# 1. Index the sorted BAM (reads mapped to the assembly, e.g. with bowtie2)
MetaCAT indexBam \${METACAT_INDEX_PARAMS} -t ${task.cpus} --bam ${bam}

# 2. Find marker-gene "seed" contigs in the assembly
MetaCAT seed \${METACAT_SEED_PARAMS} --threads ${task.cpus} \
	--fasta ${contigs} --output \${TEMP_DIR}/metacat.seed

# 3. Compute per-contig coverage (how deeply each contig is sequenced)
MetaCAT coverage \${METACAT_COVERAGE_PARAMS} --threads-count ${task.cpus} \
  --bam ${bam} --output \${TEMP_DIR}/${sample}.metacat.coverage

# 4. Bin the contigs
MetaCAT cluster \${METACAT_CLUSTER_PARAMS} --threads ${task.cpus} --fasta ${contigs} \
  --seed \${TEMP_DIR}/metacat.seed --coverage \${TEMP_DIR}/${sample}.metacat.coverage \
  --output \${TEMP_DIR}/${sample}.metacat

BIN_CONTIG_MAPPING=${sample}_bin_contig_mapping.tsv
echo -e "BIN_ID\tCONTIG\tBINNER" > \${BIN_CONTIG_MAPPING}
for bin in \$(find \${TEMP_DIR}/ -name "*fasta"); do

  ID=\$(basename \${bin} | rev | cut -d '.' -f 2 | rev)

	BIN_NAME="${sample}_bin.\${ID}.fa"

	# Append bin id to every header
	seqkit replace  -p '(.*)' -r "\\\${1} MAG=\${ID}" \$bin > \${BIN_NAME}

	# Create bin to contig mapping and add the used binner to each line
	grep ">" \${bin} | sed 's/>//g' \\
		| sed "s|^|\${BIN_NAME}\\t|g;s|\$|\\tmetacat|" >> \${BIN_CONTIG_MAPPING}
done

# return not binned fasta files
BINNED_IDS=binned.tsv
NOT_BINNED=${sample}_notBinned.fa
grep -h ">" \$(find \${TEMP_DIR}/ -name "*fasta") | tr -d ">" > \${BINNED_IDS}
if [ -s \${BINNED_IDS} ]; then
	# Get all not binned Ids
	seqkit grep -vf \${BINNED_IDS} ${contigs} \\
		| seqkit replace  -p '(.*)' -r "\\\${1} MAG=NotBinned" > \${NOT_BINNED}
else
	seqkit replace  -p '(.*)' -r "\\\${1} MAG=NotBinned" ${contigs} > \${NOT_BINNED}
fi
