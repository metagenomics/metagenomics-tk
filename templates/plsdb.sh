DATABASE=!{params.databases}/plsdb
LOCK_FILE=${DATABASE}/checksum.txt
DOWNLOAD_LINK=!{params.steps.plasmid.PLSDB.database.source}
MD5SUM=!{params.steps.plasmid.PLSDB.database.md5sum}

mkdir -p ${DATABASE}
flock ${LOCK_FILE} concurrentDownload.sh --output=${DATABASE} \
	--link=$DOWNLOAD_LINK \
	--httpsCommand="wget -O out.zip $DOWNLOAD_LINK && unzip -q out.zip && rm out.zip" \
	--s3Command="s5cmd !{params.steps?.plasmid?.PLSDB?.database?.s5cmdParams} cp ${DOWNLOAD_LINK} out.zip && unzip -q out.zip && rm out.zip" \
	--localCommand="unzip ${DOWNLOAD_LINK} -d . " \
	--expectedMD5SUM=${MD5SUM}

mash sketch -i !{params.steps.plasmid.PLSDB.additionalParams.mashSketch} -o query !{plasmids}

mash dist ${DATABASE}/out/plsdb.msh query.msh -p !{task.cpus} !{params.steps.plasmid.PLSDB.additionalParams.mashDist} > !{binID}.tsv
sort -rgk 5,5 !{binID}.tsv | sed 's/\/.*$//g' \
     | awk -v threshold=!{params.steps.plasmid.PLSDB.sharedKmerThreshold} '($5+0 > threshold) {print $0}' > matches.tsv

OUTPUT=!{binID}_kmerThreshold_!{params.steps.plasmid.PLSDB.sharedKmerThreshold}.tsv

while IFS=$"\t" read line ; do 
	hit=$(echo "$line" | cut -f 1); 
	grep -w $hit ${DATABASE}/out/plsdb.tsv \
		| sed "s/^/${line}\t/g" \
		| sed "s/^/!{sample}\t!{binID}\t/g" ; 
done <matches.tsv  > matches_detailed.tsv 

cat <(head -n 1 ${DATABASE}/out/plsdb.tsv \
	| sed "s/^/MASH_REFERENCE\tMASH_QUERY\tMASH_DISTANCE\tMASH_P_VALUE\tMASH_MATCHING_HASHES\t/g"  \
	| sed 's/^/SAMPLE\tBIN_ID\t/g') matches_detailed.tsv > ${OUTPUT}

