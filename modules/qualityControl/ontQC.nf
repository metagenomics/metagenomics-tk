nextflow.enable.dsl=2

def getOutput(SAMPLE, RUNID, TOOL, filename){
    return SAMPLE + '/' + RUNID + '/' + params.modules.qcONT.name + '/' + 
           params.modules.qcONT.version.major + "." +
           params.modules.qcONT.version.minor + "." +
           params.modules.qcONT.version.patch + 
           '/' + TOOL + '/' + filename
}

process pPorechop {

    label 'medium'

    tag "$sample"

    publishDir params.output, saveAs: { filename -> getOutput("${sample}", params.runid, "porechop", filename) }

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


process pPorechopDownload {

    label 'medium'

    tag "Sample: $sample"

    publishDir params.output, saveAs: { filename -> getOutput("${sample}", params.runid, "porechop", filename) }

    when params?.steps?.containsKey("qcONT") && params?.steps?.qcONT.containsKey("porechop")

    container "${params.porechop_image}"

    input:
    tuple val(sample), env(readUrl)

    output:
    tuple val("${sample}"), path("${sample}_qc.fq.gz"), emit: reads
    tuple file(".command.sh"), file(".command.out"), file(".command.err"), file(".command.log")

    shell:
    '''
    s5cmd !{params.steps.qcONT.porechop.download.s5cmdParams} cat ${readUrl} > reads.fq.gz 
    porechop -i reads.fq.gz -o reads.porechoped.fq.gz --threads !{task.cpus} !{params.steps.qcONT.porechop.additionalParams.porechop}
    filtlong !{params.steps.qcONT.porechop.additionalParams.filtlong} reads.porechoped.fq.gz | pigz --processes !{task.cpus} > !{sample}_qc.fq.gz
    '''
}



process pNanoPlot {

    label 'medium'

    tag "Sample: $sample"

    publishDir params.output, saveAs: { filename -> getOutput("${sample}", params.runid, "nanoplot", filename) }

    when params?.steps?.containsKey("qcONT") && params?.steps?.qcONT.containsKey("nanoplot")

    container "${params.nanoplot_image}"

    input:
    tuple val(sample), path(read1, stageAs: "reads.fq.gz")

    output:
    tuple val("${sample}"), path("*.html"), emit: html
    tuple val("${sample}"), path("*.png"), emit: png
    tuple val("${sample}"), path("*.tsv"), emit: tsv
    tuple file(".command.sh"), file(".command.out"), file(".command.err"), file(".command.log")

    shell:
    '''
    NanoPlot --threads !{task.cpus} --tsv_stats  --fastq reads.fq.gz
    csvtk -tT transpose  <(tail -n +2 NanoStats.txt) > NanoStats.tsv
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
             reads | pNanoPlot
      emit:
        reads = reads
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
            | map { it -> [ it.SAMPLE, it.READS ]} \
	    | _wONTFastq | set { results }
    emit:
      reads = results.reads
}
