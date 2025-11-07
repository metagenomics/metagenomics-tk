# if no local database is referenced, start download part
if [ -z "!{EXTRACTED_DB}" ]
then
     # polished params end with a / so no additional one is needed at the end
     DATABASE=!{params.polished.databases}mmseqs2
     mkdir -p ${DATABASE}/!{dbType}
     LOCK_FILE=${DATABASE}/!{dbType}/lock.txt

    # Check if access and secret keys are necessary for s5cmd
    if [ ! -z "!{S3_TAX_DB_ACCESS}" ]
    then
            export AWS_ACCESS_KEY_ID=!{S3_TAX_DB_ACCESS}
            export AWS_SECRET_ACCESS_KEY=!{S3_TAX_DB_SECRET}
    fi

    # Create and try to lock the given “LOCK_FILE”. If the file is locked no other process will download the same database simultaneously
    # and wait until the database is downloaded.  Depending on your database path (starts with s3:// | https:// | /../../.. )
    # one of the flock input commands is executed to download and decompress the database.
    # If a database is present at the given path, checksums are compared, if they are identical the download will be omitted.
    flock ${LOCK_FILE} concurrentDownload.sh --output=${DATABASE}/!{dbType} \
         --link=!{DOWNLOAD_LINK} \
         --httpsCommand="wgetStatic --no-check-certificate -qO- !{DOWNLOAD_LINK} | zstd -T!{task.cpus} -c -d | tar -xv " \
         --s3FileCommand="s5cmd !{S5CMD_PARAMS} cat --concurrency !{task.cpus} !{DOWNLOAD_LINK} | zstd -T!{task.cpus} -c -d | tar -xv " \
         --s3DirectoryCommand="s5cmd !{S5CMD_PARAMS} cp --concurrency !{task.cpus} !{DOWNLOAD_LINK} . " \
         --s5cmdAdditionalParams="!{S5CMD_PARAMS}" \
         --localCommand="zstd -T!{task.cpus} -c -d !{DOWNLOAD_LINK} | tar -xv " \
         --expectedMD5SUM=!{MD5SUM}

    # Path of the newly downloaded database. A mmseqs2 database consists out of multiple files,
    # the unique lookup file is searched to determine the basename of each database.
    FILEPATH=$(find ${DATABASE}/!{dbType} -name *.lookup -type f)
    # remove extension to get the files basename.
           MMSEQS2_DATABASE_DIR=${FILEPATH%.*}
else
    # If an extracted database is present use that path.
           MMSEQS2_DATABASE_DIR="!{EXTRACTED_DB}"
fi



mkdir tmp
# Fix for ownership issue https://github.com/nextflow-io/nextflow/issues/4565
chmod a+rw -R tmp

# According to the MMSeqs docs the --split-memory-limit parameter defines the RAM usage for *about* 80 percent of the total RAM consumption.
RAM_LIMIT="$(awk -v RATIO=80 -v RAM=$(echo !{task.memory} | cut -f 1 -d ' ') 'BEGIN { print int(RAM / 100 * RATIO) }')G"

# Define taxonomies
#$MMSEQS_VERSION taxonomy queryDB ${MMSEQS2_DATABASE_DIR} !{sample}_!{binType}.!{dbType}.taxresults.database tmp !{parameters}  --split-memory-limit ${RAM_LIMIT} -s !{sensitivity} --threads !{task.cpus}
# mmseqs2 searches produce output databases. These have to be converted to more useful formats.
#$MMSEQS_VERSION createtsv queryDB !{sample}_!{binType}.!{dbType}.taxresults.database !{sample}_!{binType}.!{dbType}.taxonomy.tsv --threads !{task.cpus}
#$MMSEQS_VERSION taxonomyreport ${MMSEQS2_DATABASE_DIR} !{sample}_!{binType}.!{dbType}.taxresults.database !{sample}_!{binType}.!{dbType}.krakenStyleTaxonomy.out
#$MMSEQS_VERSION taxonomyreport ${MMSEQS2_DATABASE_DIR} !{sample}_!{binType}.!{dbType}.taxresults.database !{sample}_!{binType}.!{dbType}.krona.html --report-mode 1
