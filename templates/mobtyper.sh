# Check developer documentation
MOB_TYPER_DB=""
if [ -z "!{EXTRACTED_DB}" ]
then
   DATABASE=!{params.databases}/mob_typer
   LOCK_FILE=${DATABASE}/checksum.txt
   # Download plsdb database if necessary
   mkdir -p ${DATABASE}
   flock ${LOCK_FILE} concurrentDownload.sh --output=${DATABASE} \
    --link=!{DOWNLOAD_LINK} \
    --httpsCommand=" wget -qO- !{DOWNLOAD_LINK} | tar --strip-components=1  -xvz " \
    --s3FileCommand=" s5cmd !{S5CMD_PARAMS} cat !{DOWNLOAD_LINK} | tar --strip-components=1 -xvz " \
    --s3DirectoryCommand=" s5cmd !{S5CMD_PARAMS} cp !{DOWNLOAD_LINK} mob_suite.tar.gz && tar --strip-components=1 -xvz mob_suite.tar.gz && rm mob_suite.tar.gz " \
    --s5cmdAdditionalParams="!{S5CMD_PARAMS}" \
    --localCommand="tar --strip-components=1  -xzvf !{DOWNLOAD_LINK} " \
    --expectedMD5SUM=!{MD5SUM}

   MOB_TYPER_DB=${DATABASE}/out
else
   MOB_TYPER_DB=!{EXTRACTED_DB}
fi

seqkit replace -p "\s.+" !{plasmids} \
	| seqkit seq --min-len !{MIN_LENGTH} > unzipped_plasmids.fasta

MOB_TYPER_OUT=!{sample}_!{binID}_mobtyper.tsv

# Bug fix: Mob Typer does not search for the taxa.sqlite file in the user provided database directory
sed -i " 148 i ETE3DBTAXAFILE = \"${MOB_TYPER_DB}/taxa.sqlite\"" /usr/local/lib/python3.9/dist-packages/mob_suite-3.0.3-py3.9.egg/mob_suite/constants.py

mob_typer !{ADDITIONAL_PARAMS} -d ${MOB_TYPER_DB} -n !{task.cpus} --multi --infile unzipped_plasmids.fasta --out_file out.tsv 

# Add Sample and BinID
sed  '1 s/^[^\t]*\t/CONTIG\t/' out.tsv \
  | sed -e '1 s/^/SAMPLE\tBIN_ID\t/g' -e "2,$ s/^/!{sample}\t!{binID}\t/g" > ${MOB_TYPER_OUT}
 
