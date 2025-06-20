# Create new plasmid file
seqkit range -r !{start}:!{stop} !{plasmids} > mobInput.fa

# Check developer documentation
MOB_TYPER_DB=""
if [ -z "!{EXTRACTED_DB}" ]
then
   DATABASE=!{params.databases}/mob_typer
   LOCK_FILE=${DATABASE}/checksum.txt

   # Check if access and secret keys are necessary for s5cmd
   if [ ! -z "!{S3_MobTyper_ACCESS}" ]
   then
        export AWS_ACCESS_KEY_ID=!{S3_MobTyper_ACCESS}
        export AWS_SECRET_ACCESS_KEY=!{S3_MobTyper_SECRET}
   fi

   # Download plsdb database if necessary
   mkdir -p ${DATABASE}
   flock ${LOCK_FILE} concurrentDownload.sh --output=${DATABASE} \
    --link=!{DOWNLOAD_LINK} \
    --httpsCommand=" wgetStatic --no-check-certificate -qO- !{DOWNLOAD_LINK} | tar --strip-components=1  -xvz " \
    --s3FileCommand=" s5cmd !{S5CMD_PARAMS} cat !{DOWNLOAD_LINK} | tar --strip-components=1 -xvz " \
    --s3DirectoryCommand=" s5cmd !{S5CMD_PARAMS} cp !{DOWNLOAD_LINK} . " \
    --s5cmdAdditionalParams="!{S5CMD_PARAMS}" \
    --localCommand="tar --strip-components=1  -xzvf !{DOWNLOAD_LINK} " \
    --expectedMD5SUM=!{MD5SUM}

   MOB_TYPER_DB=${DATABASE}/out
else
   MOB_TYPER_DB=!{EXTRACTED_DB}
fi

seqkit replace -p "\s.+" mobInput.fa \
	| seqkit seq --min-len !{MIN_LENGTH} > unzipped_plasmids.fasta

MOB_TYPER_OUT=!{binID}_chunk_!{start}_!{stop}_mobtyper.tsv

mob_typer !{ADDITIONAL_PARAMS} -d ${MOB_TYPER_DB} -n !{task.cpus} \
	--multi --infile unzipped_plasmids.fasta --out_file out.tsv > >(tee -a mob_typer_stdout.log) 2> >(tee -a mob_typer_stderr.log >&2)
 
# If there is not enough RAM on the host available, mobtyper exits with exit code 0 in some cases.
# The only solution is to check the stdout for exceptions.
if grep -q "Exception" mob_typer_stdout.log mob_typer_stderr.log; then
	echo "Exception found in MobTyper log.";
	exit 1 ;
else
	echo "No exception found in MobTyper log.";
fi

# Add Sample and BinID
sed  '1 s/^[^\t]*\t/CONTIG\t/' out.tsv \
  | sed -e '1 s/^/SAMPLE\tBIN_ID\t/g' -e "2,$ s/^/!{sample}\t!{binID}\t/g" > ${MOB_TYPER_OUT}
