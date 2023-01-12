nextflow.enable.dsl=2

def getOutput(SAMPLE, RUNID, TOOL, filename){
    return SAMPLE + '/' + RUNID + '/' + params.modules.qc.name + '/' + 
           params.modules.qc.version.major + "." +
           params.modules.qc.version.minor + "." +
           params.modules.qc.version.patch + 
           '/' + TOOL + '/' + filename
}


process pFastpSplit {

    label 'medium'

    tag "Sample: $sample"

    cache "deep"

    publishDir params.output, mode: "${params.publishDirMode}", saveAs: { filename -> getOutput("${sample}", params.runid, "fastp", filename) }

    when params?.steps.containsKey("qc") && params?.steps?.qc.containsKey("fastp")

    container "${params.fastp_image}"

    time params.steps.containsKey("qc") ? Utils.setTimeLimit(params.steps.qc.fastp, params.modules.qc.process.fastp.defaults, params.resources.medium) : ""

    input:
    tuple val(sample), path(read1, stageAs: "read1.fq.gz"), path(read2, stageAs: "read2.fq.gz")

    output:
    tuple val("${sample}"), path("${sample}_interleaved.qc.fq.gz"), emit: readsPair
    tuple val("${sample}"), path("${sample}_unpaired.qc.fq.gz"), emit: readsSingle
    tuple val("${sample}"), path("${sample}_fastp_summary_before.tsv"), emit: fastpSummaryBefore
    tuple val("${sample}"), path("${sample}_fastp_summary_after.tsv"), emit: fastpSummaryAfter
    tuple val("${sample}"), path("${sample}_unpaired_summary.tsv"), emit: fastpSummaryUnpaired
    tuple val("${sample}"), path("${sample}_fastp.json"), emit: fastpSummary
    tuple val("${sample}"), path("*_report.html"), emit: fastpSummaryHtml
    tuple file(".command.sh"), file(".command.out"), file(".command.err"), file(".command.log")

    shell:
    template 'fastpSplit.sh'
}


process pNonpareil {

    label 'small'

    tag "Sample: $sample"

    publishDir params.output, mode: "${params.publishDirMode}", saveAs: { filename -> getOutput("${sample}", params.runid, "nonpareil", filename) }

    when params.steps.containsKey("qc") && params?.steps?.qc.containsKey("nonpareil")

    container "${params.nonpareil_image}"

    input:
    tuple val(sample), path(interleavedReads, stageAs: 'interleaved.fq.gz'), path(unpairedReads)

    output:
    tuple file(".command.sh"), file(".command.out"), file(".command.err"), file(".command.log"), path("*.npl"), emit: log
    tuple val("${sample}"), path("*.npa"), emit: redundancyValues
    tuple val("${sample}"), path("*.npc"), emit: matesDistribution
    tuple val("${sample}"), path("*.npc"), emit: redundancySummary
    tuple val("${sample}"), path("*_nonpareil_curves.pdf"), emit: nonpareilCurves
    tuple val("${sample}"), path("*_nonpareil_index.tsv"), emit: nonpareilIndex

    shell:
    template 'nonpareil.sh'
}



process pJellyFish {

    label 'medium'

    publishDir params.output, mode: "${params.publishDirMode}", saveAs: { filename -> getOutput("${sample}", params.runid, "jellyfish", filename) }

    container "${params.jellyfish_image}"

    when params.steps.containsKey("qc") && params?.steps?.qc.containsKey("jellyfish")

    input:
    tuple val(sample), path(interleavedReads, stageAs: 'interleaved.fq.gz'), path(unpairedReads)

    output:
    tuple val("${sample}"), path("*.histo.tsv"), emit: histogram
    tuple file(".command.sh"), file(".command.out"), file(".command.err"), file(".command.log"), emit: log

    shell:
    '''
    jellyfish count !{params.steps.qc.jellyfish.additionalParams.count}   -t !{task.cpus} <(zcat !{interleavedReads} !{unpairedReads}) -o !{sample}.jf
    echo -e "FREQUENCY\tCOUNT\tSAMPLE" > !{sample}.histo.tsv
    jellyfish histo !{params.steps.qc.jellyfish.additionalParams.histo}  !{sample}.jf | sed -e 's/ /\t/g' -e "s/$/\t!{sample}/g" >> !{sample}.histo.tsv
    '''
}



process pFastpSplitDownload {

    label 'medium'

    tag "Sample: $sample"

    cache 'deep'

    publishDir params.output, mode: "${params.publishDirMode}", saveAs: { filename -> getOutput("${sample}", params.runid, "fastp", filename) }

    time params.steps.containsKey("qc") ? Utils.setTimeLimit(params.steps.qc.fastp, params.modules.qc.process.fastpDownload.defaults, params.resources.medium) : ""

    when params?.steps.containsKey("qc") && params?.steps?.qc.containsKey("fastp")

    container "${params.fastp_image}"

    input:
    tuple val(sample), env(read1Url), env(read2Url)

    output:
    tuple val("${sample}"), path("${sample}_interleaved.qc.fq.gz"), emit: readsPair
    tuple val("${sample}"), path("${sample}_unpaired.qc.fq.gz"), emit: readsSingle
    tuple val("${sample}"), path("${sample}_fastp_summary_before.tsv"), emit: fastpSummaryBefore
    tuple val("${sample}"), path("${sample}_fastp_summary_after.tsv"), emit: fastpSummaryAfter
    tuple val("${sample}"), path("${sample}_unpaired_summary.tsv"), emit: fastpSummaryUnpaired
    tuple val("${sample}"), path("${sample}_fastp.json"), emit: fastpSummary
    tuple val("${sample}"), path("*_report.html"), emit: fastpSummaryHtml
    tuple file(".command.sh"), file(".command.out"), file(".command.err"), file(".command.log")

    shell:
    template 'fastpSplitDownload.sh'
}


workflow _wFastqSplit {
       take:
         reads
       main:
            // Check if files are S3 URLs and if the download parameter is specified in the config
            FASTP_FILE_IDX = 1
            reads | branch {
              download: it[2].startsWith("s3://") && params?.steps?.qc.fastp.containsKey("download")
              noDownload: !params?.steps?.qc.fastp.containsKey("download")
             } | set { samples }

             samples.noDownload | pFastpSplit 
             samples.download | pFastpSplitDownload

             // Create summary files
             if(params.summary){
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
             }
             pFastpSplit.out.readsPair | mix(pFastpSplitDownload.out.readsPair) | set {readsPair}
             pFastpSplit.out.readsSingle | mix(pFastpSplitDownload.out.readsSingle) | set {readsSingle}

             readsPair | join(readsSingle) | (pNonpareil & pJellyFish)
      emit:
        nonpareilIndex = pNonpareil.out.nonpareilIndex
        kmerFrequencies =  pJellyFish.out.histogram 
        readsPair = readsPair
        readsSingle = readsSingle
}


/*
 * Takes a channel as input that has the following format [SAMPLE, LEFT_READ_PATH, RIGHT READ_PATH]. 
 * This module does read trimming, adapter removal and other quality control tasks.   
 * 
 * Output:
 * - Interleaved fastq files
 * - Nonpareil coverage and diversity estimation
 * - kmer frequencies
 * 
 */
workflow wShortReadQualityControlList {
  take:
    reads
  main:
    reads | _wFastqSplit | set { results } 
  emit:
    readsPair = results.readsPair
    readsSingle = results.readsSingle
    nonpareil = results.nonpareilIndex
    kmerFrequencies =  results.kmerFrequencies
}

/*
 * Takes a tab separated file of files containing reads as input and applies 
 * read trimming, adapter removal and other quality control tasks.   
 * Input file with columns seperated by tabs:
 * SAMPLE Left_Read Right_Read
 *
 * Left and right read could be https, s3 links or file path.
 * 
 * Output:
 * - Interleaved fastq files
 * - Nonpareil coverage and diversity estimation
 * - kmer frequencies
 */
workflow wShortReadQualityControlFile {
     main:
       if(!params.steps.qc.interleaved){
         Channel.from(file(params.steps.qc.input)) \
            | splitCsv(sep: '\t', header: true) \
            | map { it -> [ it.SAMPLE, it.READS1, it.READS2 ]} \
	    | _wFastqSplit | set { results }
       }
    emit:
      readsPair = results.readsPair
      readsSingle = results.readsSingle
}
