nextflow.enable.dsl=2

MODULE="assembly"
VERSION="0.1.0"
def getOutput(SAMPLE, RUNID, TOOL, filename){
    return SAMPLE + '/' + RUNID + '/' + MODULE + '/' + VERSION + '/' + TOOL + '/' + filename
}

process pFastpSplit {

    label 'large'

    tag "$sample"

    publishDir params.output, saveAs: { filename -> getOutput("${sample}", params.runid, "fastp", filename) }

    container "quay.io/biocontainers/fastp:${params.fastp_tag}"

    input:
    tuple val(sample), path(read1, stageAs: "read1.fq.gz"), path(read2, stageAs: "read2.fq.gz")

    output:
    tuple val("${sample}"), path("interleaved.fastp.fq.gz"), emit: reads_processed
    tuple val("${sample}"), path("fastp_summary_before.tsv"), emit: fastp_summary_before
    tuple val("${sample}"), path("fastp_summary_after.tsv"), emit: fastp_summary_after
    tuple val("${sample}"), path("fastp.json"), emit: fastp_summary
    tuple val("${sample}"), path("*_report.html"), emit: fastp_summary_html

    shell:
    template 'fastpSplit.sh'
}


process pMegahit {

    label 'large'

    tag "$sample"

    publishDir params.output, saveAs: { filename -> getOutput("${sample}", params.runid, "megahit", filename) }

    errorStrategy 'ignore'

    when params.steps.assembly.containsKey("megahit")

    container "vout/megahit:${params.megahit_tag}"

    input:
    tuple val(sample), path(fastqs, stageAs: 'reads.fq.gz')

    output:
    tuple val("${sample}"), path("${sample}_final.contigs.fa.gz"), emit: contigs
    tuple val("${sample}"), path("${sample}_contigs_stats.tsv"), emit: contigs_stats

    shell:
    template 'megahit.sh'
}


workflow _wFastqSplit {
       take:
         reads_table
       main:
            reads_table | splitCsv(sep: '\t', header: true) \
             | map { it -> [ it.SAMPLE, it.READS1, it.READS2 ]} \
             | pFastpSplit 

            pFastpSplit.out.reads_processed | pMegahit
             
            pFastpSplit.out.fastp_summary | collectFile(newLine: false, keepHeader: true, storeDir: params.output + "/summary/" ){ item ->
              [ "fastp_summary.tsv", item[1].text ]
            }

            pFastpSplit.out.fastp_summary_after | collectFile(newLine: false, keepHeader: true, storeDir: params.output + "/summary/" ){ item ->
              [ "fastp_summary_after.tsv", item[1].text ]
            }

            pFastpSplit.out.fastp_summary_before | collectFile(newLine: false, keepHeader: true, storeDir: params.output + "/summary/" ){ item ->
              [ "fastp_summary_before.tsv", item[1].text ]
            }

            pMegahit.out.contigs_stats | collectFile(newLine: false, keepHeader: true, storeDir: params.output + "/summary/" ){ item ->
              [ "contigs_stats.tsv", item[1].text ]
            }

      emit:
        contigs = pMegahit.out.contigs
        processed_reads = pFastpSplit.out.reads_processed
}


/*
 * Takes a tab separated file of files containing reads as input and produces assembly, 
 * binning results and contamination, completeness and taxonomy reports for mags.
 * Input file with columns seperated by tabs:
 * Dataset_ID Left_Read Right_Read
 *
 * Left and right read could be https, s3 links or file path.
 * 
 */
workflow wAssemblyFile {
     take:
       reads
     main:
       if(params.steps.assembly.interleaved){
       //     _wFastqInterleavedSeperate(reads)
       } else {
         results = _wFastqSplit(reads)
       }
    emit:
      contigs = results.contigs
      processed_reads = results.processed_reads
}
