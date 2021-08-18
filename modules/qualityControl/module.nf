nextflow.enable.dsl=2

MODULE="qc"
VERSION="0.1.0"
def getOutput(SAMPLE, RUNID, TOOL, filename){
    return SAMPLE + '/' + RUNID + '/' + MODULE + '/' + VERSION + '/' + TOOL + '/' + filename
}

process pFastpSplit {

    label 'medium'

    tag "$sample"

    publishDir params.output, saveAs: { filename -> getOutput("${sample}", params.runid, "fastp", filename) }

    when params?.steps?.qc.containsKey("fastp")

    container "quay.io/biocontainers/fastp:${params.fastp_tag}"

    input:
    tuple val(sample), path(read1, stageAs: "read1.fq.gz"), path(read2, stageAs: "read2.fq.gz")

    output:
    tuple val("${sample}"), path("${sample}_interleaved.qc.fq.gz"), emit: readsProcessed
    tuple val("${sample}"), path("fastp_summary_before.tsv"), emit: fastpSummaryBefore
    tuple val("${sample}"), path("fastp_summary_after.tsv"), emit: fastpSummaryAfter
    tuple val("${sample}"), path("fastp.json"), emit: fastpSummary
    tuple val("${sample}"), path("*_report.html"), emit: fastpSummaryHtml

    shell:
    template 'fastpSplit.sh'
}


process pFastpSplitDownload {

    label 'medium'

    echo true

    tag "$sample"

    publishDir params.output, saveAs: { filename -> getOutput("${sample}", params.runid, "fastp", filename) }

    when:
    params?.steps?.qc.containsKey("fastp")

    container "quay.io/biocontainers/fastp:${params.fastp_tag}"

    input:
    tuple val(sample), env(read1Url), env(read2Url)

    output:
    tuple val("${sample}"), path("${sample}_interleaved.qc.fq.gz"), emit: readsProcessed
    tuple val("${sample}"), path("fastp_summary_before.tsv"), emit: fastpSummaryBefore
    tuple val("${sample}"), path("fastp_summary_after.tsv"), emit: fastpSummaryAfter
    tuple val("${sample}"), path("fastp.json"), emit: fastpSummary
    tuple val("${sample}"), path("*_report.html"), emit: fastpSummaryHtml

    shell:
    template 'fastpSplitDownload.sh'
}


workflow _wFastqSplit {
       take:
         readsTable
       main:
            // Check if files are S3 URls and if the download parameter is specified config
            FASTP_FILE_IDX = 1
            readsTable | splitCsv(sep: '\t', header: true) \
             | map { it -> [ it.SAMPLE, it.READS1, it.READS2 ]} \
             | branch {
              download: it[2].startsWith("s3://") && params?.steps?.qc.containsKey("download")
              noDownload: !params?.steps?.qc.containsKey("download")
             } | set { samples }

             samples.noDownload | pFastpSplit
             samples.download | pFastpSplitDownload

             // Create summary files
             pFastpSplit.out.fastpSummary | mix(pFastpSplitDownload.out.fastpSummary) \
              | collectFile(newLine: false, keepHeader: true, storeDir: params.output + "/summary/" ){ item ->
              [ "fastp_summary.tsv", item[FASTP_FILE_IDX].text ]
             }

             pFastpSplit.out.fastpSummaryAfter | mix(pFastpSplitDownload.out.fastpSummaryAfter) \
              | collectFile(newLine: false, keepHeader: true, storeDir: params.output + "/summary/" ){ item ->
              [ "fastp_summary_after.tsv", item[FASTP_FILE_IDX].text ]
             }

             pFastpSplit.out.fastpSummaryBefore | mix(pFastpSplitDownload.out.fastpSummaryBefore) \
              | collectFile(newLine: false, keepHeader: true, storeDir: params.output + "/summary/" ){ item ->
              [ "fastp_summary_before.tsv", item[FASTP_FILE_IDX].text ]
             }
             pFastpSplit.out.readsProcessed | mix(pFastpSplitDownload.out.readsProcessed) | set {readsProcessed}
      emit:
        processed_reads = readsProcessed
}


/*
 * Takes a tab separated file of files containing reads as input and applies 
 * read trimming, adapter removal and other quality control tasks.   
 * Input file with columns seperated by tabs:
 * Dataset_ID Left_Read Right_Read
 *
 * Left and right read could be https, s3 links or file path.
 * 
 * Output:
 * Interleaved fastq files.
 * 
 */
workflow wQualityControlFile {
     take:
       readsTable
     main:
       if(params.steps.qc.interleaved){
       //     _wFastqInterleavedSeperate(reads)
       } else {
         results = _wFastqSplit(readsTable)
       }
    emit:
      processed_reads = results.processed_reads
}
