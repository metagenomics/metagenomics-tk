# if no local database is referenced, start the download part
if [ -z "!{EXTRACTED_DB}" ] 
then
# polished params end with a / so no additional one is needed at the end
    DATABASE=!{params.polished.databases}mmseqs2 
    mkdir -p ${DATABASE}/!{dbType}
    LOCK_FILE=${DATABASE}/!{dbType}/lock.txt

    # Check if access and secret keys are necessary for s5cmd
    if [ ! -z "!{S3_DB_ACCESS}" ]
    then
        export AWS_ACCESS_KEY_ID=!{S3_DB_ACCESS}
        export AWS_SECRET_ACCESS_KEY=!{S3_DB_SECRET}
    fi
         
    # Create and try to lock the given “LOCK_FILE”. If the file is locked no other process will download the same database simultaneously
    # and wait until the database is downloaded.  Depending on your database path (starts with s3:// | https:// | /../../.. ) 
    # one of the flock input commands is executed to download and decompress the database.
    # If a database is present at the given path, checksums are compared. If they are identical, the download will be omitted.  
    flock ${LOCK_FILE} concurrentDownload.sh --output=${DATABASE}/!{dbType} \
            --link=!{DOWNLOAD_LINK} \
            --httpsCommand="wget -qO- !{DOWNLOAD_LINK} | zstd -T!{task.cpus} -d -c | tar -xv " \
            --s3FileCommand="s5cmd !{S5CMD_PARAMS} cat --concurrency !{task.cpus} !{DOWNLOAD_LINK} | zstd -T!{task.cpus} -d -c | tar -xv " \
            --s3DirectoryCommand="s5cmd !{S5CMD_PARAMS} cp --concurrency !{task.cpus} !{DOWNLOAD_LINK} . " \
	    --s5cmdAdditionalParams="!{S5CMD_PARAMS}" \
            --localCommand="zstd -T!{task.cpus} -c -d !{DOWNLOAD_LINK} | tar -xv " \
            --expectedMD5SUM=!{MD5SUM}
          
    # Path of the newly downloaded database. A mmseqs2 database consists of multiple files,
    # the unique lookup file is searched to determine the basename of each database.
    FILEPATH=$(find ${DATABASE}/!{dbType} -name *.lookup -type f)
    # remove the extension to get the files basename.
    MMSEQS2_DATABASE_DIR=${FILEPATH%.*}
else
    # If an extracted database is present use that path. 
    MMSEQS2_DATABASE_DIR="!{EXTRACTED_DB}"
fi

# Check which mmseqs version is the right one
MMSEQS_VERSION=""
if grep -q "avx2" /proc/cpuinfo; then
	echo "AVX2 supported, using mmseqs_avx2"
	MMSEQS_VERSION="mmseqs_avx2"
elif grep -q "sse4_1" /proc/cpuinfo; then
	echo "SSE4.1 supported, using mmseqs_sse41"
	MMSEQS_VERSION="mmseqs_sse41"
else
	echo "Using mmseqs_sse2"
	MMSEQS_VERSION="mmseqs_sse2"
fi

MMSEQS_HEADER="query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qlen,tlen,qcov,tcov!{ADDITIONAL_COLUMNS}"
OUTPUT_TMP_TSV="!{sample}_!{binType}.!{start}.!{stop}.!{dbType}.blast.tmp.tsv"
OUTPUT_TSV="!{sample}_!{binType}.!{start}.!{stop}.!{dbType}.blast.tsv"

# According to the MMSeqs docs the --split-memory-limit parameter defines the RAM usage for *about* 80 percent of the total RAM consumption.
RAM_LIMIT="$(awk -v RATIO=80 -v RAM=$(echo !{task.memory} | cut -f 1 -d ' ') 'BEGIN { print int(RAM / 100 * RATIO) }')G"

mkdir tmp

# Create new input file
cat !{fasta} | seqkit range -r !{start}:!{stop} > input.fa

# Only mmseqs2 databases can be used for every kind of search. Inputs have to be converted first.
$MMSEQS_VERSION createdb input.fa queryDB
# Load all indices into memory to increase searching speed
$MMSEQS_VERSION touchdb --threads !{task.cpus} queryDB
$MMSEQS_VERSION touchdb --threads !{task.cpus} ${MMSEQS2_DATABASE_DIR}
$MMSEQS_VERSION search queryDB ${MMSEQS2_DATABASE_DIR} !{sample}_!{binType}.!{dbType}.results.database tmp !{parameters} \
	--threads !{task.cpus} --split-memory-limit ${RAM_LIMIT}
# mmseqs2 searches produce output databases. These have to be converted to a more useful format. The blast -outfmt 6 in this case.
$MMSEQS_VERSION convertalis queryDB ${MMSEQS2_DATABASE_DIR} !{sample}_!{binType}.!{dbType}.results.database ${OUTPUT_TMP_TSV} \
	--threads !{task.cpus} --format-output ${MMSEQS_HEADER}

# Try only if the output file is not empty, csvtk will fail otherwise
if [ -s ${OUTPUT_TMP_TSV}  ]; then
     # Add header
     sed  -i "1i $(echo ${MMSEQS_HEADER} | tr ',' '\t')" ${OUTPUT_TMP_TSV}

     # Add BIN_ID and SAMPLE column
     csvtk -t -T join -f "locus_tag;query" \
	    <(csvtk -t -T concat !{contig2GeneMapping} | csvtk -t -T cut -f SAMPLE,BIN_ID,CONTIG,locus_tag) ${OUTPUT_TMP_TSV} > ${OUTPUT_TSV}

     # use "query" column name instead of "locus_tag"
     sed -i -e "1s/locus_tag/query/" ${OUTPUT_TSV}
fi

