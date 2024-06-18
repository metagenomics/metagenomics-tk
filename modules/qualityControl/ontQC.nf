include { wSaveSettingsList } from '../config/module'

def getOutput(SAMPLE, RUNID, TOOL, filename){
    return SAMPLE + '/' + RUNID + '/' + params.modules.qcONT.name + '/' + 
           params.modules.qcONT.version.major + "." +
           params.modules.qcONT.version.minor + "." +
           params.modules.qcONT.version.patch + 
           '/' + TOOL + '/' + filename
}

/*
* This process runs Adapter trimming (Porechop) and quality control (filtlong) in one step in order to reduce disk space consumption of the work directory.
*  
*
*/
process pPorechop {

    label 'small'

    tag "$sample"

    cache 'deep'

    publishDir params.output, mode: "${params.publishDirMode}", saveAs: { filename -> getOutput("${sample}", params.runid, "porechop", filename) }

    when params?.steps?.containsKey("qcONT") && params?.steps?.qcONT.containsKey("porechop")

    container "${params.porechop_image}"

    input:
    tuple val(sample), path(read1, stageAs: "reads.fq.gz")

    output:
    tuple val("${sample}"), path("${sample}_qc.fq.gz"), emit: reads
    tuple file(".command.sh"), file(".command.out"), file(".command.err"), file(".command.log")

    shell:
    '''
    porechop -i reads.fq.gz -o reads.porechoped.fq.gz --threads !{task.cpus} !{params.steps.qcONT.porechop.additionalParams.porechop}
    filtlong !{params.steps.qcONT.porechop.additionalParams.filtlong} reads.porechoped.fq.gz | pigz --processes !{task.cpus} > !{sample}_qc.fq.gz
    '''
}


/*
* See Porechop
*/
process pPorechopDownload {

    label 'small'

    tag "Sample: $sample"

    cache 'deep'

    publishDir params.output, mode: "${params.publishDirMode}", saveAs: { filename -> getOutput("${sample}", params.runid, "porechop", filename) }

    when params?.steps?.containsKey("qcONT") && params?.steps?.qcONT.containsKey("porechop")

    container "${params.porechop_image}"

    input:
    tuple val(sample), env(readUrl)

    output:
    tuple val("${sample}"), path("${sample}_qc.fq.gz"), emit: reads
    tuple file(".command.sh"), file(".command.out"), file(".command.err"), file(".command.log")

    shell:
    '''
    set -o pipefail
    
    s5cmd !{params.steps.qcONT.porechop.download.s5cmdParams} cat ${readUrl} > reads.fq.gz 2> error.log

    # This if statement solves issue https://github.com/pbelmann/meta-omics-toolkit/issues/166
    if grep -q "reset by peer" error.log; then
       echo "Network issue found. Exiting with exit code 1";
       exit 1 ;
    else
       echo "No network issue found";
    fi

    porechop -i reads.fq.gz -o reads.porechoped.fq.gz --threads !{task.cpus} !{params.steps.qcONT.porechop.additionalParams.porechop}
    filtlong !{params.steps.qcONT.porechop.additionalParams.filtlong} reads.porechoped.fq.gz \
	| pigz --processes !{task.cpus} > !{sample}_qc.fq.gz
    '''
}


/*
*
* This process removes sequences that originate from human DNA.
*
*/
process pFilterHumanONT {

    label 'medium'

    tag "Sample: $sample"

    publishDir params.output, mode: "${params.publishDirMode}", saveAs: { filename -> getOutput("${sample}", params.runid, "filterHumanONT", filename) }

    when params.steps.containsKey("qcONT") && params?.steps?.qc.containsKey("filterHumanONT")

    container "${params.scrubber_image}"

    input:
    tuple val(sample), path(longReads)

    output:
    tuple file(".command.sh"), file(".command.out"), file(".command.err"), file(".command.log"), emit: log
    tuple val("${sample}"), path("*_filtered.fq.gz"), emit: filteredSeqs
    tuple val("${sample}"), path("*_summary_before.tsv"), emit: summaryBefore
    tuple val("${sample}"), path("*_summary_after.tsv"), emit: summaryAfter
    tuple val("${sample}"), path("*_removed.fq.gz"), optional: true, emit: removed

    shell:
    EXTRACTED_DB=params.steps?.qcONT?.filterHumanONT?.database?.extractedDBPath ?: ""
    DOWNLOAD_LINK=params.steps?.qcONT?.filterHumanONT?.database?.download?.source ?: ""
    MD5SUM=params?.steps?.qcONT?.filterHumanONT?.database?.download?.md5sum ?: ""
    S5CMD_PARAMS=params.steps?.qcONT?.filterHumanONT?.database?.download?.s5cmd?.params ?: ""
    output = getOutput("${sample}", params.runid, "filterHumanONT", "")
    ADDITIONAL_PARAMS=params.steps?.qcONT?.filterHumanONT?.additionalParams ?: ""
    S3_filter_ACCESS=params?.steps?.qcONT?.filterHumanONT?.database?.download?.s5cmd && S5CMD_PARAMS.indexOf("--no-sign-request") == -1 ? "\$S3_filter_ACCESS" : ""
    S3_filter_SECRET=params?.steps?.qcONT?.filterHumanONT?.database?.download?.s5cmd && S5CMD_PARAMS.indexOf("--no-sign-request") == -1 ? "\$S3_filter_SECRET" : ""
    template 'filterHumanONT.sh'
}


process pNanoPlot {

    label 'small'

    tag "Sample: $sample"

    publishDir params.output, mode: "${params.publishDirMode}", saveAs: { filename -> getOutput("${sample}", params.runid, "nanoplot", filename) }

    when params?.steps?.containsKey("qcONT") && params?.steps?.qcONT.containsKey("nanoplot")

    container "${params.nanoplot_image}"

    input:
    tuple val(sample), path(read1, stageAs: "reads.fq.gz")

    output:
    tuple val("${sample}"), path("*.html"), emit: html
    tuple val("${sample}"), path("*.png"), emit: png
    tuple val("${sample}"), path("*.tsv"), emit: tsv
    tuple val("${sample}"), env(MEDIAN_QUALITY), emit: medianQuality
    tuple file(".command.sh"), file(".command.out"), file(".command.err"), file(".command.log")

    shell:
    '''
    NanoPlot --threads !{task.cpus} --tsv_stats  --fastq reads.fq.gz
    csvtk -tT transpose  <(tail -n +2 NanoStats.txt)  > NanoStatsTmp.tsv
    paste -d$'\t' <(echo -e "SAMPLE\n!{sample}") NanoStatsTmp.tsv > NanoStats.tsv
    MEDIAN_QUALITY=$(cut -f 8 NanoStats.tsv | tail -n 1)
    '''
}



workflow _wONTFastq {
       take:
         reads
       main:
            // Check if files are S3 URLs and if the download parameter is specified in the config
            FASTP_FILE_IDX = 1
            reads | branch {
              download: it[FASTP_FILE_IDX].startsWith("s3://") && params?.steps?.qcONT.porechop.containsKey("download")
              noDownload: !params?.steps?.qcONT.porechop.containsKey("download")
             } | set { samples }

             samples.noDownload | pPorechop
             samples.download | pPorechopDownload

             pPorechop.out.reads | mix(pPorechopDownload.out.reads)  | set { reads }

             filteredSeqs = Channel.empty()
             if(params.steps.containsKey("qcONT") && params.steps.qcONT.containsKey("filterHumanONT")){
             	reads | pFilterHumanONT | set{ filteredSeqs }
	     } else {
                reads | set { filteredSeqs }
             }

             filteredSeqs | pNanoPlot
      emit:
        reads = reads
        medianQuality = pNanoPlot.out.medianQuality
}


/*
 * Takes a channel as input that has the following format [SAMPLE,  READ_PATH]. 
 * This module does read trimming, adapter removal and other quality control tasks.   
 * 
 * Output:
 * Interleaved fastq files.
 * 
 */
workflow wOntQualityControlList {
  take:
    reads
  main:
    reads | _wONTFastq | set { results } 
  emit:
    reads = results.reads
    medianQuality = results.medianQuality
}

/*
 * Takes a tab separated file of files containing reads as input and applies 
 * read trimming, adapter removal and other quality control tasks.   
 * Input file with columns seperated by tabs:
 * SAMPLE Reads
 *
 * Reads could be provided as https, s3 links or file path.
 * 
 * Output:
 * Quality controlled reads.
 * 
 */
workflow wOntQualityControlFile {
     take:
       readsTable
     main:
        Channel.from(file(params.steps.qcONT.input)) 
            | splitCsv(sep: '\t', header: true) \
            | map { it -> [ it.SAMPLE, it.READS ]} | set { reads }

	_wONTFastq(reads) | set { results }

        SAMPLE_IDX = 0
        wSaveSettingsList(reads | map { it -> it[SAMPLE_IDX] })
    emit:
      reads = results.reads
}
