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

    containerOptions Utils.getDockerMount(params?.steps?.plasmid?.ViralVerifyPlasmid?.database, params)

    container "${params.viralVerify_image}"

    beforeScript "mkdir -p ${params.polished.databases}"

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
    '''
    FILTER_STRING="!{params.steps?.plasmid?.ViralVerifyPlasmid?.filterString}"

    pigz -p !{task.cpus} -fdc !{plasmids} > unzipped.fasta

    PFAM_FILE=""
    if [ -z "!{EXTRACTED_DB}" ]
    then
         DATABASE=!{params.polished.databases}/pfam
         LOCK_FILE=${DATABASE}/lock.txt

         mkdir -p ${DATABASE}
         flock ${LOCK_FILE} concurrentDownload.sh --output=${DATABASE} \
            --link=!{DOWNLOAD_LINK} \
            --httpsCommand="wget -qO- !{DOWNLOAD_LINK} | pigz -fdc > pfam.hmm " \
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

    beforeScript "mkdir -p ${params.polished.databases}"

    containerOptions Utils.getDockerMount(params.steps?.plasmid?.MobTyper?.database, params)

    input:
    tuple val(sample), val(binID), path(plasmids), val(start), val(stop)

    output:
    tuple val("${binID}_chunk_${start}_${stop}"), val("${output}"), val(params.LOG_LEVELS.INFO), file(".command.sh"), \
      file(".command.out"), file(".command.err"), file(".command.log"), emit: logs
    tuple val("${sample}"), val("${binID}"), val("MobTyper"), \
	path("${binID}_chunk_${start}_${stop}_mobtyper.tsv"), emit: plasmidsStats, optional: true

    shell:
    EXTRACTED_DB=params.steps?.plasmid?.MobTyper?.database?.extractedDBPath ?: ""
    DOWNLOAD_LINK=params?.steps?.plasmid?.MobTyper?.database?.download?.source ?: ""
    MD5SUM=params?.steps?.plasmid?.MobTyper?.database?.download?.md5sum ?: ""
    S5CMD_PARAMS=params.steps?.plasmid?.MobTyper?.database?.download?.s5cmd?.params ?: ""
    ADDITIONAL_PARAMS=params.steps?.plasmid?.MobTyper?.additionalParams ?: ""
    MIN_LENGTH=params.steps?.plasmid?.MobTyper?.minLength ?: ""
    output = getOutput("${sample}", params.runid, "MobTyper", "")
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
    tuple val(sample), val(binID), path(assembly), val(start), val(stop)

    output:
    tuple val("${binID}_chunk_${start}_${stop}"), val("${output}"), val(params.LOG_LEVELS.INFO), file(".command.sh"), \
      file(".command.out"), file(".command.err"), file(".command.log"), emit: logs
    tuple val("${sample}"), val("${binID}"), val("PlasClass"), \
	path("${binID}_chunk_${start}_${stop}_plasClass.tsv"), emit: probabilities, optional: true

    shell:
    output = getOutput("${sample}", params.runid, "PlasClass", "")
    filter = params?.steps?.plasmid?.PlasClass?.threshold ? "TRUE" : ""
    template("plasClass.sh")
}


process pPlaton {

    label 'highmemMedium'

    tag "Sample: $sample, BinId: $binID"

    containerOptions " --user root:root " + Utils.getDockerMount(params.steps?.plasmid?.Platon?.database, params)

    publishDir params.output, mode: "${params.publishDirMode}", saveAs: { filename -> getOutput("${sample}", params.runid, "Platon", filename) }, \
        pattern: "{**.tsv}"

    when params.steps.containsKey("plasmid") && params.steps.plasmid?.containsKey("Platon")

    container "${params.platon_image}"

    beforeScript "mkdir -p ${params.polished.databases}"

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
    '''
    # --meta is not supported by Platon
    sed -i "458 i       '-p', 'meta', " /usr/local/lib/python3.9/site-packages/platon/functions.py

    # Check developer documentation
    PLATON_DB=""
    if [ -z "!{EXTRACTED_DB}" ]
    then 
      DATABASE=!{params.polished.databases}/platon
      LOCK_FILE=${DATABASE}/checksum.txt

      # Download plsdb database if necessary
      mkdir -p ${DATABASE}
      flock ${LOCK_FILE} concurrentDownload.sh --output=${DATABASE} \
	--link=!{DOWNLOAD_LINK} \
	--httpsCommand="wget -qO- !{DOWNLOAD_LINK} | tar -xvz " \
	--s3FileCommand="s5cmd !{S5CMD_PARAMS} cat !{DOWNLOAD_LINK} | tar -xvz  " \
	--s3DirectoryCommand="s5cmd !{S5CMD_PARAMS} cp --concurrency !{task.cpus} !{DOWNLOAD_LINK} . " \
	--s5cmdAdditionalParams="!{S5CMD_PARAMS}" \
	--localCommand="tar xzvf !{DOWNLOAD_LINK} " \
	--expectedMD5SUM=!{MD5SUM}

      PLATON_DB=${DATABASE}/out/db
    else
      PLATON_DB=!{EXTRACTED_DB}
    fi

    pigz -p !{task.cpus} -fdc !{assembly} > assembly.fasta

    # In some cases prodigal fails because of a too short query sequence. In such cases the process should end with exit code 0.
    trap 'if [ "$?" == 1 ] && ( grep -q "ORFs failed" assembly.log || grep -q "ORFs=0$" assembly.log || grep -q "Error detecting input file format. First line seems to be blank." assembly.log ); then echo "Protein Prediction Failed"; exit 0; fi' EXIT
    platon assembly.fasta !{ADDITIONAL_PARAMS} --db ${PLATON_DB} --mode sensitivity -t !{task.cpus}

    if [ -n "$(find . -name '*.tsv')" ]; then
      # Add Sample and BinId
      sed -e '1 s/^/SAMPLE\tBIN_ID\t/g' -e "1 s/\tID\t/\tCONTIG\t/"  -e "2,$ s/^/!{sample}\t!{binID}\t/g" *.tsv > !{binID}_platon.tsv
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
       for file in !{contigHeaderFiles}; do 
         TMP_FILE=tmp/${file}
         FINAL_FILE=header/${file}
         MISSING_FILE=missing/${file}
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



