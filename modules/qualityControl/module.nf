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
    tuple val("${sample}"), path("${sample}_interleaved.qc.fq.gz"), emit: reads_processed
    tuple val("${sample}"), path("fastp_summary_before.tsv"), emit: fastp_summary_before
    tuple val("${sample}"), path("fastp_summary_after.tsv"), emit: fastp_summary_after
    tuple val("${sample}"), path("fastp.json"), emit: fastp_summary
    tuple val("${sample}"), path("*_report.html"), emit: fastp_summary_html

    shell:
    template 'fastpSplit.sh'
}



workflow _wFastqSplit {
       take:
         readsTable
       main:
            FASTP_FILE_IDX = 1
            readsTable | splitCsv(sep: '\t', header: true) \
             | map { it -> [ it.SAMPLE, it.READS1, it.READS2 ]} \
             | pFastpSplit

            pFastpSplit.out.fastp_summary | collectFile(newLine: false, keepHeader: true, storeDir: params.output + "/summary/" ){ item ->
              [ "fastp_summary.tsv", item[FASTP_FILE_IDX].text ]
            }

            pFastpSplit.out.fastp_summary_after | collectFile(newLine: false, keepHeader: true, storeDir: params.output + "/summary/" ){ item ->
              [ "fastp_summary_after.tsv", item[FASTP_FILE_IDX].text ]
            }

            pFastpSplit.out.fastp_summary_before | collectFile(newLine: false, keepHeader: true, storeDir: params.output + "/summary/" ){ item ->
              [ "fastp_summary_before.tsv", item[FASTP_FILE_IDX].text ]
            }
      emit:
        processed_reads = pFastpSplit.out.reads_processed
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
