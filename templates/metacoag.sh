
# Activate Conda environment
conda init bash
source /opt/conda/bin/activate
conda activate metacoag

# compute contig abundance 
coverm contig -t !{task.cpus} --bam-files !{bam} | tail -n +2 > assembly_depth.tsv

# Update contig names to make it consistent with the gfa and flyes assembly info file
csvtk replace -r {kv}  -f 1 -H -t -k <(awk -F $'\t' ' { t = $1; $1 = $2; $2 = t; print; }' OFS=$'\t' !{headerMapping} ) assembly_depth.tsv -p "(.*)"  > assembly_depth_final.tsv
seqkit replace -p '(.*)' -r '{kv}' \
	-k <(awk -F $'\t' ' { t = $1; $1 = $2; $2 = t; print; } ' OFS=$'\t' !{headerMapping})  !{contigs} > contigs.fa.gz

# Run MetaCoAG
metacoag --assembler flye --graph !{graph}  --nthreads !{task.cpus} --contigs contigs.fa.gz \
	--paths !{flyeAssemblyInfo} !{params.steps.binning.metacoag.additionalParams} --abundance assembly_depth_final.tsv --output .

# Rename contig and file names according to toolkit specification
BIN_CONTIG_MAPPING=!{sample}_bin_contig_mapping.tsv
BINNED_IDS=binned.tsv
for bin in  $(find bins -name "bin*.fasta") ; do

	# Get id of the bin (e.g get 2 of the bin SAMPLEID_bin.2.fa)
	ID=$(echo ${bin} | cut -f 2 -d '_' | cut -f 1 -d '.')

	BIN_NAME="!{sample}_bin.${ID}.fa"

        # Add BIN ID
	seqkit replace -p '(.*)'  -r '{kv}' -k !{headerMapping} ${bin} \
	| tee >( grep -h ">" | tr -d '>' >> ${BINNED_IDS} ) \
 	| seqkit replace  -p '(.*)' -r "\${1} MAG=${ID}" > ${BIN_NAME}

	# Create bin to contig mapping
	grep ">" ${bin} | sed 's/>//g' \
        | csvtk replace -r {kv}  -f 1 -H -t -k !{headerMapping} -p "(.*)" \
     	| sed "s/^/${BIN_NAME}\t/g" >> ${BIN_CONTIG_MAPPING}
done

# return not binned fasta files
BINNED_IDS=binned.tsv
NOT_BINNED=!{sample}_notBinned.fa
grep -h ">" $(basename !{contigs})*/bin* | tr -d ">" > ${BINNED_IDS}
if [ -s ${BINNED_IDS} ]; then
	# Get all not binned Ids
	seqkit grep -vf ${BINNED_IDS} !{contigs} \
	| seqkit replace  -p '(.*)' -r "\${1} MAG=NotBinned" > ${NOT_BINNED}
else
	seqkit replace  -p '(.*)' -r "\${1} MAG=NotBinned" !{contigs} > ${NOT_BINNED}
fi

