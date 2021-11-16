DATABASE=!{params.databases}/plsdb
CHECKSUM_FILE=${DATABASE}/checksum.txt
DOWNLOAD_LINK=!{params.steps.plasmid.PLSDB.database.source}
MD5SUM=!{params.steps.plasmid.PLSDB.database.md5sum}

mkdir -p ${DATABASE}
flock ${CHECKSUM_FILE} concurrentDownload.sh --output=${DATABASE} \
	--checkpoint=${CHECKSUM_FILE} \
	--command="wget -O out.zip $DOWNLOAD_LINK && unzip -q out.zip && rm out.zip" \
	--mode=MD5SUM --expectedMD5SUM=${MD5SUM}

mash sketch -i !{params.steps.plasmid.PLSDB.additionalParams.mashSketch} -o query !{plasmids}

mash dist ${DATABASE}/out/plsdb.msh query.msh -p !{task.cpus} !{params.steps.plasmid.PLSDB.additionalParams.mashDist} > !{binID}.tsv
sort -rgk 5,5 !{binID}.tsv | sed 's/\/.*$//g' \
     | awk -v threshold=!{params.steps.plasmid.PLSDB.sharedKmerThreshold} '($5+0 > threshold) {print $0}' | cut -f 1  > ids.tsv

OUTPUT=!{binID}_kmerThreshold_!{params.steps.plasmid.PLSDB.sharedKmerThreshold}.tsv

cat  <(head -n 1 ${DATABASE}/out/plsdb.tsv | sed 's/^/SAMPLE\tBIN_ID\t/g') \
	<(grep -f ids.tsv ${DATABASE}/out/plsdb.tsv | sed "s/^/!{sample}\t!{binID}\t/g" ) > ${OUTPUT} 
