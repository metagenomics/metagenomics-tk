nextflow.enable.dsl=2

include { pDumpLogs } from '../utils/processes'

def getOutput(SAMPLE, RUNID, TOOL, filename){
    return SAMPLE + '/' + RUNID + '/' + params.modules.plasmids.name + '/' +
          params.modules.plasmids.version.major + "." +
          params.modules.plasmids.version.minor + "." +
          params.modules.plasmids.version.patch +
          '/' + TOOL + '/' + filename
}


process pViralVerifyPlasmid {

    label 'small'

    tag "Sample: $sample, BinID: $binID"

    publishDir params.output, mode: "${params.publishDirMode}", saveAs: { filename -> getOutput("${sample}", params.runid, "ViralVerifyPlasmid", filename) }, \
        pattern: "{**.tsv}"

    when params.steps.containsKey("plasmid") && params.steps.plasmid?.containsKey("ViralVerifyPlasmid")

    shell = ['/bin/bash']

    containerOptions Utils.getDockerMount(params?.steps?.plasmid?.ViralVerifyPlasmid?.database, params) + Utils.getDockerNetwork()

    container "${params.viralVerify_image}"

    beforeScript Utils.getCreateDatabaseDirCommand("${params.polished.databases}")

    secret { "${S3_ViralVerifyPlasmid_ACCESS}"!="" ? ["S3_ViralVerifyPlasmid_ACCESS", "S3_ViralVerifyPlasmid_SECRET"] : [] } 

    input:
    tuple val(sample), val(binID), path(plasmids)

    output:
    tuple val("${sample}"), val("${binID}"), val("ViralVerifyPlasmid"), path("${binID}_viralverifyplasmid.tsv"), emit: plasmidsStats, optional: true
    tuple val("${binID}"), val("${output}"), val(params.LOG_LEVELS.INFO), file(".command.sh"), \
      file(".command.out"), file(".command.err"), file(".command.log"), emit: logs

    shell:
    EXTRACTED_DB=params.steps?.plasmid?.ViralVerifyPlasmid?.database?.extractedDBPath ?: ""
    DOWNLOAD_LINK=params.steps?.plasmid?.ViralVerifyPlasmid?.database?.download?.source ?: ""
    MD5SUM=params?.steps?.plasmid?.ViralVerifyPlasmid?.database?.download?.md5sum ?: ""
    S5CMD_PARAMS=params.steps?.plasmid?.ViralVerifyPlasmid?.database?.download?.s5cmd?.params ?: ""
    output = getOutput("${sample}", params.runid, "ViralVerifyPlasmid", "")
    ADDITIONAL_PARAMS=params.steps?.plasmid?.ViralVerifyPlasmid?.additionalParams ?: ""
    S3_ViralVerifyPlasmid_ACCESS=params?.steps?.plasmid?.ViralVerifyPlasmid?.database?.download?.s5cmd && S5CMD_PARAMS.indexOf("--no-sign-request") == -1 ? "\$S3_ViralVerifyPlasmid_ACCESS" : ""
    S3_ViralVerifyPlasmid_SECRET=params?.steps?.plasmid?.ViralVerifyPlasmid?.database?.download?.s5cmd && S5CMD_PARAMS.indexOf("--no-sign-request") == -1 ? "\$S3_ViralVerifyPlasmid_SECRET" : ""
    '''
    FILTER_STRING="!{params.steps?.plasmid?.ViralVerifyPlasmid?.filterString}"

    pigz -p !{task.cpus} -fdc !{plasmids} > unzipped.fasta

    PFAM_FILE=""
    if [ -z "!{EXTRACTED_DB}" ]
    then
         DATABASE=!{params.polished.databases}/pfam
         LOCK_FILE=${DATABASE}/lock.txt

	 # Check if access and secret keys are necessary for s5cmd
         if [ ! -z "!{S3_ViralVerifyPlasmid_ACCESS}" ]
         then
             export AWS_ACCESS_KEY_ID=!{S3_ViralVerifyPlasmid_ACCESS}
             export AWS_SECRET_ACCESS_KEY=!{S3_ViralVerifyPlasmid_SECRET}
         fi

         mkdir -p ${DATABASE}
         flock ${LOCK_FILE} concurrentDownload.sh --output=${DATABASE} \
            --link=!{DOWNLOAD_LINK} \
            --httpsCommand="wgetStatic --no-check-certificate -qO- !{DOWNLOAD_LINK} | pigz -fdc > pfam.hmm " \
            --s3FileCommand="s5cmd !{S5CMD_PARAMS} cat --concurrency !{task.cpus} !{DOWNLOAD_LINK} | pigz -fdc > pfam.hmm  " \
            --s3DirectoryCommand="s5cmd !{S5CMD_PARAMS} cp --concurrency !{task.cpus} !{DOWNLOAD_LINK} . && mv * pfam.hmm " \
            --s5cmdAdditionalParams="!{S5CMD_PARAMS}" \
            --localCommand="gunzip -c !{DOWNLOAD_LINK} > ./pfam.hmm " \
            --expectedMD5SUM=!{MD5SUM}

          PFAM_FILE="${DATABASE}/out/pfam.hmm"
    else
          PFAM_FILE="!{EXTRACTED_DB}"
    fi

    # ViralVerify should exit gracefully if the query sequence is too short for hmmsearch
    trap 'if [[ $? == 1 && ! -s unzipped*proteins_circ.fa ]]; then echo "hmmsearch could not detect any protein"; exit 0; fi' EXIT
    viralverify !{ADDITIONAL_PARAMS}  --hmm ${PFAM_FILE} -p -t !{task.cpus} -f unzipped.fasta -o .

    TMP_OUTPUT="!{binID}_viralverifyplasmid_tmp.tsv"
    FINAL_OUTPUT="!{binID}_viralverifyplasmid.tsv"

    # In case that no sequences could be found then at least an empty header file should be created
    echo -e "CONTIG\tSAMPLE\tBIN_ID" > ${FINAL_OUTPUT} 

    if [ -n "$(find . -name '*.csv')" ]; then
      # Filter viral verify output by user provided strings 
      # and add Sample and BinID
      sed  's/,/\t/g' *.csv \
	| sed -E -n "1p;/${FILTER_STRING}/p" \
        | sed -e '1 s/^[^\t]*\t/CONTIG\t/' -e '1 s/^/SAMPLE\tBIN_ID\t/g' -e "2,$ s/^/!{sample}\t!{binID}\t/g" > ${TMP_OUTPUT}
      if [ $(wc -l < ${TMP_OUTPUT}) -gt 1 ]; then
        mv ${TMP_OUTPUT} ${FINAL_OUTPUT}
      fi
    fi
    '''
}


process pMobTyper {

    label 'highmemMedium'

    tag "Sample: $sample, BinID: $binID"

    publishDir params.output, mode: "${params.publishDirMode}", saveAs: { filename -> getOutput("${sample}", params.runid, "MobTyper", filename) }, \
        pattern: "{**.tsv}"

    when params.steps.containsKey("plasmid") && params.steps.plasmid?.containsKey("MobTyper")

    container "${params.mobSuite_image}"

    beforeScript Utils.getCreateDatabaseDirCommand("${params.polished.databases}")

    containerOptions Utils.getDockerMount(params.steps?.plasmid?.MobTyper?.database, params) + Utils.getDockerNetwork()

    secret { "${S3_MobTyper_ACCESS}"!="" ? ["S3_MobTyper_ACCESS", "S3_MobTyper_SECRET"] : [] } 

    input:
    tuple val(sample), val(binID), path(plasmids), val(start), val(stop), val(chunkSize)

    output:
    tuple val("${binID}_chunk_${start}_${stop}"), val("${output}"), val(params.LOG_LEVELS.INFO), file(".command.sh"), \
      file(".command.out"), file(".command.err"), file(".command.log"), emit: logs
    tuple val("${sample}"), val("${binID}"), val("MobTyper"), \
	path("${binID}_chunk_${start}_${stop}_mobtyper.tsv"), val("${chunkSize}"), emit: plasmidsStats, optional: true

    shell:
    EXTRACTED_DB=params.steps?.plasmid?.MobTyper?.database?.extractedDBPath ?: ""
    DOWNLOAD_LINK=params?.steps?.plasmid?.MobTyper?.database?.download?.source ?: ""
    MD5SUM=params?.steps?.plasmid?.MobTyper?.database?.download?.md5sum ?: ""
    S5CMD_PARAMS=params.steps?.plasmid?.MobTyper?.database?.download?.s5cmd?.params ?: ""
    ADDITIONAL_PARAMS=params.steps?.plasmid?.MobTyper?.additionalParams ?: ""
    MIN_LENGTH=params.steps?.plasmid?.MobTyper?.minLength ?: ""
    output = getOutput("${sample}", params.runid, "MobTyper", "")
    S3_MobTyper_ACCESS=params.steps?.plasmid?.MobTyper?.database?.download?.s5cmd && S5CMD_PARAMS.indexOf("--no-sign-request") == -1 ? "\$S3_MobTyper_ACCESS" : ""
    S3_MobTyper_SECRET=params.steps?.plasmid?.MobTyper?.database?.download?.s5cmd && S5CMD_PARAMS.indexOf("--no-sign-request") == -1 ? "\$S3_MobTyper_SECRET" : ""
    template("mobtyper.sh")
}


process pPlasClass {

    label 'highmemMedium'

    tag "Sample: $sample, BinID: $binID"

    publishDir params.output, mode: "${params.publishDirMode}", saveAs: { filename -> getOutput("${sample}", params.runid, "PlasClass", filename) }, \
        pattern: "{**.tsv}"

    when params.steps.containsKey("plasmid") && params.steps.plasmid?.containsKey("PlasClass")

    container "${params.PlasClass_image}"

    input:
    tuple val(sample), val(binID), path(assembly), val(start), val(stop), val(chunkSize)

    output:
    tuple val("${binID}_chunk_${start}_${stop}"), val("${output}"), val(params.LOG_LEVELS.INFO), file(".command.sh"), \
      file(".command.out"), file(".command.err"), file(".command.log"), emit: logs
    tuple val("${sample}"), val("${binID}"), val("PlasClass"), \
	path("${binID}_chunk_${start}_${stop}_plasClass.tsv"), val("${chunkSize}"), emit: probabilities, optional: true

    shell:
    output = getOutput("${sample}", params.runid, "PlasClass", "")
    filter = params?.steps?.plasmid?.PlasClass?.threshold ? "TRUE" : ""
    template("plasClass.sh")
}


process pPlaton {

    label 'highmemMedium'

    tag "Sample: $sample, BinId: $binID"

    containerOptions " --user root:root " + Utils.getDockerMount(params.steps?.plasmid?.Platon?.database, params) + Utils.getDockerNetwork()

    publishDir params.output, mode: "${params.publishDirMode}", saveAs: { filename -> getOutput("${sample}", params.runid, "Platon", filename) }, \
        pattern: "{**.tsv}"

    when params.steps.containsKey("plasmid") && params.steps.plasmid?.containsKey("Platon")

    container "${params.platon_image}"

    beforeScript Utils.getCreateDatabaseDirCommand("${params.polished.databases}")

    secret { "${S3_Platon_ACCESS}"!="" ? ["S3_Platon_ACCESS", "S3_Platon_SECRET"] : [] } 

    input:
    tuple val(sample), val(binID), path(assembly)

    output:
    tuple val("${sample}"), val("${binID}"), val("Platon"), path("${binID}_platon.tsv"), optional: true, emit: plasmidsStats
    tuple val("${binID}"), val("${output}"), val(params.LOG_LEVELS.INFO), file(".command.sh"), \
      file(".command.out"), file(".command.err"), file(".command.log"), emit: logs

    shell:
    output = getOutput("${sample}", params.runid, "Platon", "")
    ADDITIONAL_PARAMS=params.steps?.plasmid?.Platon?.additionalParams ?: ""
    EXTRACTED_DB=params.steps?.plasmid?.Platon?.database?.extractedDBPath ?: ""
    DOWNLOAD_LINK=params?.steps?.plasmid?.Platon?.database?.download?.source ?: ""
    MD5SUM=params?.steps?.plasmid?.Platon?.database?.download?.md5sum ?: ""
    S5CMD_PARAMS=params.steps?.plasmid?.Platon?.database?.download?.s5cmd?.params ?: ""
    S3_Platon_ACCESS=params.steps?.plasmid?.Platon?.database?.download?.s5cmd && S5CMD_PARAMS.indexOf("--no-sign-request") == -1 ? "\$S3_Platon_ACCESS" : ""
    S3_Platon_SECRET=params.steps?.plasmid?.Platon?.database?.download?.s5cmd && S5CMD_PARAMS.indexOf("--no-sign-request") == -1 ? "\$S3_Platon_SECRET" : ""
    '''
    # --meta is not supported by Platon
    sed -i "458 i       '-p', 'meta', " /usr/local/lib/python3.9/site-packages/platon/functions.py

    # Check developer documentation
    PLATON_DB=""
    if [ -z "!{EXTRACTED_DB}" ]
    then 
      DATABASE=!{params.polished.databases}/platon
      LOCK_FILE=${DATABASE}/checksum.txt

      # Check if access and secret keys are necessary for s5cmd
      if [ ! -z "!{S3_Platon_ACCESS}" ]
      then
          export AWS_ACCESS_KEY_ID=!{S3_Platon_ACCESS}
          export AWS_SECRET_ACCESS_KEY=!{S3_Platon_SECRET}
      fi

      # Download plsdb database if necessary
      mkdir -p ${DATABASE}
      flock ${LOCK_FILE} concurrentDownload.sh --output=${DATABASE} \
	--link=!{DOWNLOAD_LINK} \
	--httpsCommand="wgetStatic --no-check-certificate -qO- !{DOWNLOAD_LINK} | tar -xvz " \
	--s3FileCommand="s5cmd !{S5CMD_PARAMS} cat !{DOWNLOAD_LINK} | tar -xvz  " \
	--s3DirectoryCommand="mkdir db && cd db && s5cmd !{S5CMD_PARAMS} cp --concurrency !{task.cpus} !{DOWNLOAD_LINK} . " \
	--s5cmdAdditionalParams="!{S5CMD_PARAMS}" \
	--localCommand="tar xzvf !{DOWNLOAD_LINK} " \
	--expectedMD5SUM=!{MD5SUM}

      PLATON_DB=${DATABASE}/out/db
    else
      PLATON_DB=!{EXTRACTED_DB}
    fi

    FINAL_OUTPUT=!{binID}_platon.tsv

    pigz -p !{task.cpus} -fdc !{assembly} > assembly.fasta

    # In some cases prodigal fails because of a too short query sequence. In such cases the process should end with exit code 0.
    trap 'if [ "$?" == 1 ] && ( grep -q "ORFs failed" assembly.log || grep -q "ORFs=0$" assembly.log || grep -q "Error detecting input file format. First line seems to be blank." assembly.log ); then echo "Protein Prediction Failed"; exit 0; fi' EXIT
    platon assembly.fasta !{ADDITIONAL_PARAMS} --db ${PLATON_DB} --mode sensitivity -t !{task.cpus}

    if [ -n "$(find . -name '*.tsv')" ]; then
      # Add Sample and BinId
      sed -e '1 s/^/SAMPLE\tBIN_ID\t/g' -e "1 s/\tID\t/\tCONTIG\t/"  -e "2,$ s/^/!{sample}\t!{binID}\t/g" *.tsv > ${FINAL_OUTPUT}
    else
      # In case that no sequences could be found then at least an empty header file should be created
      echo -e "SAMPLE\tBIN_ID\tCONTIG" > ${FINAL_OUTPUT} 
    fi
    '''
}


process pFilter {

    label 'small'

    tag "Sample: $sample, BinId: $binID"

    container "${params.ubuntu_image}"

    publishDir params.output, mode: "${params.publishDirMode}", saveAs: { filename -> getOutput("${sample}", params.runid, "filtered", filename) }, \
        pattern: "{**.tsv,**.fasta.gz}"

    when params.steps.containsKey("plasmid") && params.steps.plasmid.containsKey("Filter")

    input:
    tuple val(sample), val(binID), val(size), path(contigs), path(contigHeaderFiles)

    output:
    tuple val("${sample}"), val("${binID}"), path("${binID}_filtered.fasta.gz"), emit: plasmids, optional: true
    tuple val("${sample}"), val("${binID}"), path("${binID}_filtered.tsv"), emit: plasmidsTsv, optional: true
    tuple val("${sample}"), val("${binID}"), path("${binID}_detection_tools.tsv"), emit: detectionTools, optional: true
    tuple val("${binID}"), val("${output}"), val(params.LOG_LEVELS.INFO), file(".command.sh"), \
        file(".command.out"), file(".command.err"), file(".command.log"), emit: logs

    shell:
    MIN_LENGTH=params.steps?.plasmid?.Filter?.minLength
    NUMBER_OF_CONTIGS=size+1
    switch(params.steps?.plasmid?.Filter.method) {
      case "OR":
       '''
       mkdir header
       mkdir tmp
       mkdir missing
       mkdir input

       for file in !{contigHeaderFiles}; do
	  if [ $(wc -l < ${file}) -gt 1 ]; then
            mv ${file} input/${file}
          fi
       done

       for file in $(ls -1 input/*); do 
         BASENAME=$(basename ${file})
         TMP_FILE=tmp/${BASENAME}
         FINAL_FILE=header/${BASENAME}
         MISSING_FILE=missing/${BASENAME}
         METHOD="$(echo ${file} | rev | cut -d '_' -f 1 | rev | cut -d '.' -f 1)"
         csvtk cut -f CONTIG --tabs ${file} | sed -e "2,$ s/$/\tTRUE/g"  -e "1 s/$/\t${METHOD}/g" > ${TMP_FILE}
	 csvtk cut -f CONTIG --tabs ${file} | tail -n +2 >> filtered_header.tsv
         if [ $(wc -l < ${TMP_FILE}) -gt 1 ]; then
                mv ${TMP_FILE} ${FINAL_FILE}
         else
                mv ${TMP_FILE} ${MISSING_FILE}
         fi
       done

       PLASMID_OUT_FASTA=!{binID}_filtered.fasta.gz 
       PLASMID_OUT_TSV=!{binID}_filtered.tsv

       if [ -s filtered_header.tsv ]; then
         sort filtered_header.tsv | uniq > filtered_sorted_header.tsv

         csvtk -t join -f 1 <(echo "CONTIG"; seqkit fx2tab --name --only-id !{contigs}) header/*  -k --na FALSE > !{binID}_detection_tools.tsv

         for file in $(ls -1 missing/) ; do  
            MISS_METHOD=$(csvtk cut -f -CONTIG --tabs $file);  
            sed -i -e "2,$ s/$/\tFALSE/g" -e "1 s/$/\t${MISS_METHOD}/g" !{binID}_detection_tools.tsv
         done

         seqkit grep -f filtered_sorted_header.tsv !{contigs} | seqkit seq --min-len !{MIN_LENGTH} \
         | pigz -c > ${PLASMID_OUT_FASTA}

         seqkit fx2tab -H --length --only-id --gc --name ${PLASMID_OUT_FASTA} > ${PLASMID_OUT_TSV}
       fi
       '''
       break;
      case "AND":
       template("filterAnd.sh")
       break;
    }
}



