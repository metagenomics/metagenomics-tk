nextflow.enable.dsl=2

include { pDumpLogs } from '../utils/processes'

def getOutput(SAMPLE, RUNID, TOOL, filename){
    return SAMPLE + '/' + RUNID + '/' + params.modules.plasmids.name + '/' +
          params.modules.plasmids.version.major + "." +
          params.modules.plasmids.version.minor + "." +
          params.modules.plasmids.version.patch +
          '/' + TOOL + '/' + filename
}


process pSCAPP {

    label 'medium'

    tag "Sample: $sample, BinID: $binID"

    publishDir params.output, saveAs: { filename -> getOutput("${sample}", params.runid, "SCAPP", filename) }

    when params.steps.containsKey("plasmid") && params.steps.plasmid?.containsKey("SCAPP")

    container "${params.SCAPP_image}"

    input:
    tuple val(sample), path(assemblyGraph), val(maxKmer), path(bam)

    output:
    tuple val("${sample}"), path("${sample}_plasmids.fasta.gz"), emit: plasmids, optional: true
    tuple file(".command.sh"), file(".command.out"), file(".command.err"), file(".command.log")
    tuple val("${sample}"), path("${sample}_plasmids_stats.tsv"), emit: plasmidsStats, optional: true

    shell:
    template("scapp.sh")
}


process pViralVerifyPlasmid {

    label 'medium'

    tag "Sample: $sample, BinID: $binID"

    publishDir params.output, saveAs: { filename -> getOutput("${sample}", params.runid, "ViralVerifyPlasmid", filename) }, \
        pattern: "{**.tsv}"

    when params.steps.containsKey("plasmid") && params.steps.plasmid?.containsKey("ViralVerifyPlasmid")

    containerOptions Utils.getDockerMount(params?.steps?.plasmid?.ViralVerifyPlasmid?.database, params)

    container "${params.viralVerify_image}"

    input:
    tuple val(sample), val(binID), path(plasmids)

    output:
    tuple val("${sample}"), val("${binID}"), val("ViralVerifyPlasmid"), path("${sample}_${binID}_viralverifyplasmid.tsv"), emit: plasmidsStats, optional: true
    tuple val("${sample}_${binID}"), val("${output}"), val(params.LOG_LEVELS.INFO), file(".command.sh"), \
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
            --httpsCommand="wget -O - !{DOWNLOAD_LINK} | pigz -fdc > pfam.hmm " \
            --s3FileCommand="s5cmd !{S5CMD_PARAMS} cat !{DOWNLOAD_LINK} | pigz -fdc > pfam.hmm  " \
            --s3DirectoryCommand="s5cmd !{S5CMD_PARAMS} cp !{DOWNLOAD_LINK} pfam.hmm.gz && pigz -fdc pfam.hmm.gz > pfam.hmm && rm pfam.hmm.gz " \
            --s5cmdAdditionalParams="!{S5CMD_PARAMS}" \
            --localCommand="gunzip -c !{DOWNLOAD_LINK} > ./pfam.hmm " \
            --expectedMD5SUM=!{MD5SUM}

          PFAM_FILE="${DATABASE}/out/pfam.hmm"
    else
          PFAM_FILE="!{EXTRACTED_DB}"
    fi

    viralverify !{ADDITIONAL_PARAMS}  --hmm ${PFAM_FILE} -p -t !{task.cpus} -f unzipped.fasta -o .

    # Filter viral verify output by user provided strings 
    # and add Sample and BinID
    sed  's/,/\t/g' *.csv \
	| sed -E -n "1p;/${FILTER_STRING}/p" \
        | sed -e '1 s/^[^\t]*\t/CONTIG\t/' -e '1 s/^/SAMPLE\tBIN_ID\t/g' -e "2,$ s/^/!{sample}\t!{binID}\t/g"  > !{sample}_!{binID}_viralverifyplasmid.tsv
    '''
}


process pMobTyper {

    label 'large'

    tag "Sample: $sample, BinID: $binID"

    publishDir params.output, saveAs: { filename -> getOutput("${sample}", params.runid, "MobTyper", filename) }, \
        pattern: "{**.tsv}"

    when params.steps.containsKey("plasmid") && params.steps.plasmid?.containsKey("MobTyper")

    container "${params.mobSuite_image}"

    containerOptions Utils.getDockerMount(params.steps?.plasmid?.MobTyper?.database, params)

    input:
    tuple val(sample), val(binID), path(plasmids)

    output:
    tuple val("${sample}_${binID}"), val("${output}"), val(params.LOG_LEVELS.INFO), file(".command.sh"), \
      file(".command.out"), file(".command.err"), file(".command.log"), emit: logs
    tuple val("${sample}"), val("${binID}"), val("MobTyper"), \
	path("${sample}_${binID}_mobtyper.tsv"), emit: plasmidsStats, optional: true

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

    label 'medium'

    tag "$sample $binID"

    publishDir params.output, saveAs: { filename -> getOutput("${sample}", params.runid, "PlasClass", filename) }, \
        pattern: "{**.tsv}"

    when params.steps.containsKey("plasmid") && params.steps.plasmid?.containsKey("PlasClass")

    container "${params.PlasClass_image}"

    input:
    tuple val(sample), val(binID), path(assembly)

    output:
    tuple val("${sample}"), val("${binID}"), val("PlasClass"), path("${sample}_${binID}_plasclass.tsv"), emit: probabilities
    tuple val("${sample}_${binID}"), val("${output}"), val(params.LOG_LEVELS.INFO), file(".command.sh"), \
        file(".command.out"), file(".command.err"), file(".command.log"), emit: logs

    shell:
    output = getOutput("${sample}", params.runid, "PlasClass", "")
    filter = params?.steps?.plasmid?.PlasClass?.threshold ? "TRUE" : ""
    template("plasClass.sh")
}


process pPlaton {

    label 'medium'

    tag "Sample: $sample, BinId: $binID"

    containerOptions " --user root:root " + Utils.getDockerMount(params.steps?.plasmid?.Platon?.database, params)

    publishDir params.output, saveAs: { filename -> getOutput("${sample}", params.runid, "Platon", filename) }, \
        pattern: "{**.tsv}"

    when params.steps.containsKey("plasmid") && params.steps.plasmid?.containsKey("Platon")

    container "${params.platon_image}"

    input:
    tuple val(sample), val(binID), path(assembly)

    output:
    tuple val("${sample}"), val("${binID}"), val("Platon"), path("${sample}_${binID}_platon.tsv"), emit: plasmidsStats
    tuple val("${sample}_${binID}"), val("${output}"), val(params.LOG_LEVELS.INFO), file(".command.sh"), \
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
      DATABASE=!{params.databases}/platon
      LOCK_FILE=${DATABASE}/checksum.txt

      # Download plsdb database if necessary
      mkdir -p ${DATABASE}
      flock ${LOCK_FILE} concurrentDownload.sh --output=${DATABASE} \
	--link=!{DOWNLOAD_LINK} \
	--httpsCommand="wget -qO- !{DOWNLOAD_LINK} | tar -xvz " \
	--s3FileCommand="s5cmd !{S5CMD_PARAMS} cat !{DOWNLOAD_LINK} | tar -xvz  " \
	--s3DirectoryCommand="s5cmd !{S5CMD_PARAMS} cp !{DOWNLOAD_LINK} platon.tar.gz && tar -xzvf platon.tar.gz && rm platon.tar.gz " \
	--s5cmdAdditionalParams="!{S5CMD_PARAMS}" \
	--localCommand="tar xzvf !{DOWNLOAD_LINK} " \
	--expectedMD5SUM=!{MD5SUM}

      PLATON_DB=${DATABASE}/out/db
    else
      PLATON_DB=!{EXTRACTED_DB}
    fi

    pigz -p !{task.cpus} -fdc !{assembly} > assembly.fasta
    platon assembly.fasta !{ADDITIONAL_PARAMS} --db ${PLATON_DB} --mode sensitivity -t !{task.cpus}

    # Add Sample and BinId
    sed -e '1 s/^/SAMPLE\tBIN_ID\t/g' -e "2,$ s/^/!{sample}\t!{binID}\t/g" *.tsv > !{sample}_!{binID}_platon.tsv
    '''
}


process pFilter {

    label 'small'

    tag "$sample $binID"

    container "${params.ubuntu_image}"

    publishDir params.output, saveAs: { filename -> getOutput("${sample}", params.runid, "filtered", filename) }, \
        pattern: "{**.tsv,**.fasta.gz}"

    when params.steps.containsKey("plasmid") && params.steps.plasmid.containsKey("Filter")

    input:
    tuple val(sample), val(binID), val(size), path(contigs), path(contigHeaderFiles)

    output:
    tuple val("${sample}"), val("${binID}"), path("${sample}_${binID}_plasmids_filtered.fasta.gz"), emit: plasmids, optional: true
    tuple val("${sample}_${binID}"), val("${output}"), val(params.LOG_LEVELS.INFO), file(".command.sh"), \
        file(".command.out"), file(".command.err"), file(".command.log"), emit: logs

    shell:
    MIN_LENGTH=params.steps?.plasmid?.Filter?.minLength
    NUMBER_OF_CONTIGS=size+1
    switch(params.steps?.plasmid?.Filter.method) {
      case "OR":
       '''
       for file in !{contigHeaderFiles}; do 
    	 csvtk cut -f CONTIG --tabs ${file} | tail -n +2 >> filtered_header.tsv
       done
       sort filtered_header.tsv | uniq > filtered_sorted_header.tsv
       seqkit grep -f filtered_sorted_header.tsv !{contigs} | seqkit seq --min-len !{MIN_LENGTH} \
	| pigz -c > !{sample}_!{binID}_plasmids_filtered.fasta.gz
       '''
       break;
      case "AND":
       template("filter_and.sh")
       break;
    }
}



