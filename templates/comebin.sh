if [[ "$COMP_TYPE" == "cuda" ]]; then
		# Check developer documentation
		if [ -z "${EXTRACTED_DB}" ]
		then
			DATABASE=${params.databases}/checkm
			LOCK_FILE=\${DATABASE}/lock.txt

			export CHECKM_DATA_PATH=${params.databases}/checkm/out
			echo '{"dataRoot": "!{params.databases}/checkm/out", "remoteManifestURL": "https://data.ace.uq.edu.au/public/CheckM_databases/", "manifestType": "CheckM", "remoteManifestName": ".dmanifest", "localManifestName": ".dmanifest"}' > /tmp/DATA_CONFIG

			if [ ! -z "${S3_checkm_ACCESS}" ]
			then
				export AWS_ACCESS_KEY_ID=${S3_checkm_ACCESS}
				export AWS_SECRET_ACCESS_KEY=${S3_checkm_SECRET}
			fi

			# Download checkm database if necessary
			mkdir -p \${DATABASE}
			flock \${LOCK_FILE} concurrentDownload.sh --output=\${DATABASE} \
				--link=${DOWNLOAD_LINK} \
				--httpsCommand="wgetStatic --no-check-certificate -qO- ${DOWNLOAD_LINK} | tar -xzv " \
				--s3FileCommand="s5cmd ${S5CMD_PARAMS} cat --concurrency ${task.cpus} ${DOWNLOAD_LINK} | tar -xzv " \
				--s3DirectoryCommand="s5cmd ${S5CMD_PARAMS} cp --concurrency ${task.cpus} ${DOWNLOAD_LINK} . " \
				--s5cmdAdditionalParams="${S5CMD_PARAMS}" \
				--localCommand="tar -xzvf ${DOWNLOAD_LINK}" \
				--expectedMD5SUM=${MD5SUM}
		else
			export CHECKM_DATA_PATH=${EXTRACTED_DB}
			echo '{"dataRoot": "!{EXTRACTED_DB}", "remoteManifestURL": "https://data.ace.uq.edu.au/public/CheckM_databases/", "manifestType": "CheckM", "remoteManifestName": ".dmanifest", "localManifestName": ".dmanifest"}' > /tmp/DATA_CONFIG
		fi
	
		COMP_SCRIPT=/comebin_env/.pixi/envs/default/bin/run_comebin.sh
else
		COMP_SCRIPT=/usr/local/bin/run_comebin.sh
fi


COMEBIN_PARAMS="${params.steps.binning.comebin.additionalParams.comebin}"
MIN_LENGTH="${params.steps.binning.comebin.additionalParams.length}"

mkdir old_bins

zcat ${contigs}  \
	| seqkit seq --min-len \${MIN_LENGTH} > contigs_unzipped.fa

sed -i \
  -e 's/realpath -e --/[ -e "\${OPTARG}" ] \\&\\& realpath/g' \
  -e 's/realpath -m --/mkdir -p "\$(dirname "\${OPTARG}")" \\&\\& realpath/g' \
  -e 's/realpath --/realpath/g' \
  \${COMP_SCRIPT}

run_comebin.sh -a contigs_unzipped.fa \
-d ${COMP_TYPE} \
-o old_bins \
-p bam \
-t ${task.cpus} \
\${COMEBIN_PARAMS}

BIN_CONTIG_MAPPING=${sample}_bin_contig_mapping.tsv
echo -e "BIN_ID\tCONTIG\tBINNER" > \${BIN_CONTIG_MAPPING}
for bin in \$(find old_bins/comebin_res/comebin_res_bins -name "*.fa"); do
	BIN_NAME="${sample}_bin.\$(basename \${bin})"

	# Get id of the bin (e.g get 2 of the bin SAMPLEID_bin.2.fa)
	ID=\$(echo \${BIN_NAME} | rev | cut -d '.' -f 2 | rev)

	# Append bin id to every header
	seqkit replace  -p '(.*)' -r "\\\${1} MAG=\${ID}" \$bin > \${BIN_NAME}

	# Create bin to contig mapping and add the used binner to each line
	grep ">" \${bin} | sed 's/>//g' \\
		| sed "s|^|\${BIN_NAME}\\t|g;s|\$|\\tcomebin|" >> \${BIN_CONTIG_MAPPING}
done

# return not binned fasta files
BINNED_IDS=binned.tsv
NOT_BINNED=${sample}_notBinned.fa
grep -h ">" \$(find old_bins -name "bin*.fa") | tr -d ">" > \${BINNED_IDS}
if [ -s \${BINNED_IDS} ]; then
	# Get all not binned Ids
	seqkit grep -vf \${BINNED_IDS} ${contigs} \\
		| seqkit replace  -p '(.*)' -r "\\\${1} MAG=NotBinned" > \${NOT_BINNED}
else
	seqkit replace  -p '(.*)' -r "\\\${1} MAG=NotBinned" ${contigs} > \${NOT_BINNED}
fi
