# if no local database is referenced, start the download part
if [ -z "${EXTRACTED_DB}" ]
then
    # polished params end with a / so no additional one is needed at the end
    DATABASE=${params.polished.databases}diamond
    mkdir -p \${DATABASE}/${dbType}
    LOCK_FILE=\${DATABASE}/${dbType}/lock.txt

    # Check if access and secret keys are necessary for s5cmd
    if [ ! -z "${S3_DB_ACCESS}" ]
    then
        export AWS_ACCESS_KEY_ID=${S3_DB_ACCESS}
        export AWS_SECRET_ACCESS_KEY=${S3_DB_SECRET}
    fi

    # DIAMOND databases are distributed as gzipped .dmnd files (see filterHuman). The files are decompressed into the database directory,
    # so the checksum is computed over the decompressed database. Files that are not compressed are used as they are.
    # gzip of the DIAMOND image (busybox) can not pass uncompressed files through, therefore the file is tested first.
    DECOMPRESS="if gzip -t db.download 2> /dev/null; then gzip -dc db.download > db.dmnd; rm db.download; else mv db.download db.dmnd; fi"
    flock \${LOCK_FILE} concurrentDownload.sh --output=\${DATABASE}/${dbType} \
            --link=${DOWNLOAD_LINK} \
            --httpsCommand="wgetStatic --no-check-certificate -qO db.download ${DOWNLOAD_LINK} ; \${DECOMPRESS}" \
            --s3FileCommand="s5cmd ${S5CMD_PARAMS} cp --concurrency ${task.cpus} ${DOWNLOAD_LINK} db.download ; \${DECOMPRESS}" \
            --s3DirectoryCommand="s5cmd ${S5CMD_PARAMS} cp --concurrency ${task.cpus} ${DOWNLOAD_LINK} . && ( ls *.gz > /dev/null 2>&1 && gunzip -f *.gz || true )" \
            --s5cmdAdditionalParams="${S5CMD_PARAMS}" \
            --localCommand="cp ${DOWNLOAD_LINK} db.download ; \${DECOMPRESS}" \
            --expectedMD5SUM=${MD5SUM}

    # A DIAMOND database consists of a single .dmnd file.
    DIAMOND_DATABASE=\$(find \${DATABASE}/${dbType}/out -name "*.dmnd" -type f | head -n 1)
else
    # If an extracted database is present use that path.
    DIAMOND_DATABASE="${EXTRACTED_DB}"
fi

DIAMOND_HEADER="query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qlen,tlen,qcov,tcov${ADDITIONAL_COLUMNS}"
DIAMOND_FORMAT="qseqid,sseqid,pident,length,mismatch,gapopen,qstart,qend,sstart,send,evalue,bitscore,qlen,slen,qcovhsp,scovhsp${ADDITIONAL_COLUMNS}"
OUTPUT_TMP_TSV="${sample}_${binType}.${start}.${stop}.${dbType}.blast.tmp.tsv"
OUTPUT_TSV="${sample}_${binType}.${start}.${stop}.${dbType}.blast.tsv"

mkdir tmp

# Create new input file
cat ${fasta} | seqkit range -r ${start}:${stop} > input.fa

diamond blastp --db \${DIAMOND_DATABASE} --query input.fa --out \${OUTPUT_TMP_TSV} \
	--threads ${task.cpus} --tmpdir tmp --outfmt 6 \$(echo \${DIAMOND_FORMAT} | tr ',' ' ') ${parameters}

# Try only if the output file is not empty, csvtk will fail otherwise
if [ -s \${OUTPUT_TMP_TSV} ]; then
     # MMseqs2 reports coverages as fractions (0-1), DIAMOND as percentages (0-100). The columns are converted
     # so that the output of both tools is comparable (qcov and tcov are the columns 15 and 16).
     awk 'BEGIN { FS=OFS="\\t" } { \$15 = \$15 / 100; \$16 = \$16 / 100; print }' \${OUTPUT_TMP_TSV} > \${OUTPUT_TMP_TSV}.converted
     mv \${OUTPUT_TMP_TSV}.converted \${OUTPUT_TMP_TSV}

     # Add header
     sed -i "1i \$(echo \${DIAMOND_HEADER} | tr ',' '\\t')" \${OUTPUT_TMP_TSV}

     # Add BIN_ID and SAMPLE column
     csvtk -t -T join -f "locus_tag;query" \
	    <(csvtk -t -T concat ${contig2GeneMapping} | csvtk -t -T cut -f SAMPLE,BIN_ID,CONTIG,locus_tag) \${OUTPUT_TMP_TSV} > \${OUTPUT_TSV}

     # use "query" column name instead of "locus_tag"
     sed -i -e "1s/locus_tag/query/" \${OUTPUT_TSV}
fi
