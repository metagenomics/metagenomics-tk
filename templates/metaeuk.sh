if [ -z "!{EXTRACTED_DB}" ] 
then
# polished params end with a / so no additional one is needed at the end
    DATABASE=!{params.polished.databases}metaeuk 
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
            --httpsCommand="wgetStatic --no-check-certificate -qO- !{DOWNLOAD_LINK} | zstd -T!{task.cpus} -d -c | tar -xv " \
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

OUTPUT_TMP_GFF="!{sample}_!{binType}.!{start}.!{stop}.!{dbType}.tmp.gff"
OUTPUT_GFF="!{sample}_!{binType}.!{start}.!{stop}.!{dbType}.gff"
OUTPUT_TMP_FAA="!{sample}_!{binType}.!{start}.!{stop}.!{dbType}.tmp.faa"
OUTPUT_FAA="!{sample}_!{binType}.!{start}.!{stop}.!{dbType}.faa"
OUTPUT_TMP_TAX_PER_PRED="!{sample}_!{binType}.!{start}.!{stop}.!{dbType}.tmp.tax_per_pred.tsv"
OUTPUT_TAX_PER_PRED="!{sample}_!{binType}.!{start}.!{stop}.!{dbType}.tax_per_pred.tsv"
OUTPUT_TMP_TAX_PER_CONTIG="!{sample}_!{binType}.!{start}.!{stop}.!{dbType}.tmp.tax_per_contig.tsv"
OUTPUT_TAX_PER_CONTIG="!{sample}_!{binType}.!{start}.!{stop}.!{dbType}.tax_per_contig.tsv"
OUTPUT_TMP_CODON="!{sample}_!{binType}.!{start}.!{stop}.!{dbType}.tmp.codon.fas"
OUTPUT_CODON="!{sample}_!{binType}.!{start}.!{stop}.!{dbType}.codon.fas"

# According to the MMSeqs docs the --split-memory-limit parameter defines the RAM usage for *about* 80 percent of the total RAM consumption.
RAM_LIMIT="$(awk -v RATIO=80 -v RAM=$(echo !{task.memory} | cut -f 1 -d ' ') 'BEGIN { print int(RAM / 100 * RATIO) }')G"

mkdir tmp
# Fix for ownership issue https://github.com/nextflow-io/nextflow/issues/4565
chmod a+rw -R tmp

# Create new input file
cat !{fasta} | seqkit range -r !{start}:!{stop} > input.fa

# Check which mmseqs version is the right one
MMSEQS_VERSION=""
if grep -q "avx2" /proc/cpuinfo; then
	echo "AVX2 supported, using metaeuk_avx2"
	MMSEQS_VERSION="metaeuk_avx2"
elif grep -q "sse4_1" /proc/cpuinfo; then
	echo "SSE4.1 supported, using metaeuk_sse41"
	MMSEQS_VERSION="metaeuk_sse41"
else
	echo "Using metaeuk_sse2"
	MMSEQS_VERSION="metaeuk_sse2"
fi

${MMSEQS_VERSION} easy-predict input.fa ${MMSEQS2_DATABASE_DIR} predsResults tempFolder !{parameters} \
    --threads !{task.cpus} --split-memory-limit ${RAM_LIMIT}

${MMSEQS_VERSION} taxtocontig tempFolder/latest/contigs  predsResults.fas predsResults.headersMap.tsv \
    ${MMSEQS2_DATABASE_DIR} taxResults tmpTaxFolder --majority 0.5 --tax-lineage 1 --lca-mode 2

mv predsResults.gff ${OUTPUT_TMP_GFF}
mv predsResults.fas ${OUTPUT_TMP_FAA}
mv taxResults_tax_per_pred.tsv  ${OUTPUT_TMP_TAX_PER_PRED}
mv taxResults_tax_per_contig.tsv  ${OUTPUT_TMP_TAX_PER_CONTIG}
mv predsResults.codon.fas ${OUTPUT_TMP_CODON}

if [ -s ${OUTPUT_TMP_GFF}  ]; then

     GFF_HEADER="seqname,source,feature,start,end,score,strand,frame,attribute"
     # Add header
     sed  -i "1i $(echo ${GFF_HEADER} | tr ',' '\t')" ${OUTPUT_TMP_GFF}

     # Add BIN_ID column
     csvtk -t -T join -f "CONTIG;seqname" \
	    <(csvtk -t -T concat !{contig2GeneMapping} | csvtk -t -T cut -f BIN_ID,CONTIG) ${OUTPUT_TMP_GFF} > ${OUTPUT_GFF}
fi

if [ -s ${OUTPUT_TMP_FAA}  ]; then
     # Add BIN_ID column
     mkdir faaTMP
     seqkit fx2tab ${OUTPUT_TMP_FAA} | cut -d$'\t' -f 1 > faaTMP/header.tsv
     paste -d$'\t' faaTMP/header.tsv  <(cut -d '|' -f 2 faaTMP/header.tsv) | sed '1i HEADER\tCONTIG' > faaTMP/header_contig.tsv

     csvtk -t -T join -f "CONTIG;CONTIG" \
         <(csvtk -t -T concat !{contig2GeneMapping} | csvtk -t -T cut -f BIN_ID,CONTIG) \
         faaTMP/header_contig.tsv > faaTMP/header_contig_sample.tsv 
     
     awk -v OFS='\t' '{print $3, $3, $2, $1}' faaTMP/header_contig_sample.tsv > faaTMP/header_contig_sample_reordered.tsv
     sed -i -e 's/\t/|||||/g' -e 's/|||||/\t/1' faaTMP/header_contig_sample_reordered.tsv 
     seqkit replace -p '^(\S+)' -r '{kv}' -k faaTMP/header_contig_sample_reordered.tsv ${OUTPUT_TMP_FAA} > ${OUTPUT_FAA}
fi


if [ -s ${OUTPUT_TMP_CODON}  ]; then
     # Add BIN_ID column
     #
     mkdir fasTMP
     seqkit fx2tab ${OUTPUT_TMP_CODON}| cut -d$'\t' -f 1 > fasTMP/header.tsv
     paste -d$'\t' fasTMP/header.tsv  <(cut -d '|' -f 2 fasTMP/header.tsv) | sed '1i HEADER\tCONTIG' > fasTMP/header_contig.tsv

     csvtk -t -T join -f "CONTIG;CONTIG" \
         <(csvtk -t -T concat !{contig2GeneMapping} | csvtk -t -T cut -f BIN_ID,CONTIG) \
         fasTMP/header_contig.tsv > fasTMP/header_contig_sample.tsv 
     
     awk -v OFS='\t' '{print $3, $3, $2, $1}' fasTMP/header_contig_sample.tsv > fasTMP/header_contig_sample_reordered.tsv
     sed -i -e 's/\t/|||||/g' -e 's/|||||/\t/1' fasTMP/header_contig_sample_reordered.tsv 
     seqkit replace -p '^(\S+)' -r '{kv}' -k fasTMP/header_contig_sample_reordered.tsv ${OUTPUT_TMP_CODON} > ${OUTPUT_CODON}
fi

# Try only if the output file is not empty, csvtk will fail otherwise
if [ -s ${OUTPUT_TMP_TAX_PER_PRED} ]; then
     TAX_PER_PRED_HEADER="QUERY,TAX_ID,RANK_NAME,RANK_SCIENTIFIC_NAME,LINEAGE"

     mkdir taxPredTMP

     # Add header
     sed  -i "1i $(echo ${TAX_PER_PRED_HEADER} | tr ',' '\t')" ${OUTPUT_TMP_TAX_PER_PRED}

     paste -d '\t' <(csvtk -t -T cut -f QUERY ${OUTPUT_TMP_TAX_PER_PRED} | cut -d '|' -f 2 | sed '1 s/QUERY/CONTIG/') \
         ${OUTPUT_TMP_TAX_PER_PRED} > taxPredTMP/taxResults_tax_per_pred_output.tsv

     csvtk -t -T join -f "CONTIG" \
	    <(csvtk -t -T concat !{contig2GeneMapping} \
	    | csvtk -t -T cut -f BIN_ID,CONTIG) taxPredTMP/taxResults_tax_per_pred_output.tsv >  ${OUTPUT_TAX_PER_PRED}
fi

if [ -s ${OUTPUT_TMP_TAX_PER_CONTIG} ]; then
     TAX_PER_CONTIG_HEADER="CONTIG,TAX_ID,RANK_NAME,RANK_SCIENTIFIC_NAME,RETAINED,TAXONOMICALLY_ASSIGNED,AGREE_WITH_CONTIG_LABEL,SUPPORT_RECEIVED,LINEAGE"

     sed  -i "1i $(echo ${TAX_PER_CONTIG_HEADER} | tr ',' '\t')" ${OUTPUT_TMP_TAX_PER_CONTIG}

     csvtk -t -T join -f "CONTIG" \
	    <(csvtk -t -T concat !{contig2GeneMapping} \
	    | csvtk -t -T cut -f BIN_ID,CONTIG) ${OUTPUT_TMP_TAX_PER_CONTIG} > ${OUTPUT_TAX_PER_CONTIG}
fi
