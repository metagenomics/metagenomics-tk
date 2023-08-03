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




process pKMC {

    label 'small'

    tag "Sample: $sample"

    publishDir params.output, mode: "${params.publishDirMode}", saveAs: { filename -> getOutput("${sample}", params.runid, "jellyfish", filename) }

    container "${params.kmc_image}"

    time params.steps.containsKey("qc") && params?.steps?.qc.containsKey("kmc") ? Utils.setTimeLimit(params.steps.qc.kmc, params.modules.qc.process.kmc.defaults, params.resources.small) : ""

    when params.steps.containsKey("qc") && params?.steps?.qc.containsKey("kmc")

    input:
    tuple val(sample), path(interleavedReads, stageAs: 'interleaved.fq.gz'), path(unpairedReads)

    output:
    tuple val("${sample}"), path("*.71.histo.tsv"), path("*.21.histo.tsv"), path("*.13.histo.tsv"), emit: histogram
    tuple val("${sample}"), path("*.json"), emit: details
    tuple file(".command.sh"), file(".command.out"), file(".command.err"), file(".command.log"), emit: log

    shell:
    '''
    mkdir work
    cat !{interleavedReads} !{unpairedReads} > input.fq.gz

    kmc -j!{sample}.13.kmc.json !{params.steps.qc.kmc.additionalParams.count} -m$(echo !{task.memory} | cut -d ' ' -f 1) -t!{task.cpus} -ci2 -k13  input.fq.gz 13mers work
    kmc_tools -t!{task.cpus} transform 13mers histogram !{sample}.13.histo.tmp.tsv !{params.steps.qc.kmc.additionalParams.histo}

    rm -rf 13mers

    kmc -j!{sample}.21.kmc.json !{params.steps.qc.kmc.additionalParams.count} -m$(echo !{task.memory} | cut -d ' ' -f 1) -t!{task.cpus} -ci2 -k21 input.fq.gz 21mers work
    kmc_tools -t!{task.cpus} transform 21mers histogram !{sample}.21.histo.tmp.tsv !{params.steps.qc.kmc.additionalParams.histo}

    rm -rf 21mers

    kmc -j!{sample}.71.kmc.json !{params.steps.qc.kmc.additionalParams.count} -m$(echo !{task.memory} | cut -d ' ' -f 1) -t!{task.cpus} -ci2 -k71 input.fq.gz 71mers work
    kmc_tools -t!{task.cpus} transform 71mers histogram !{sample}.71.histo.tmp.tsv !{params.steps.qc.kmc.additionalParams.histo}

    rm -rf 71mers

    echo -e "FREQUENCY\tCOUNT\tSAMPLE" > !{sample}.21.histo.tsv
    cat !{sample}.21.json | jq -r '.Stats."#k-mers_below_min_threshold"' \
	| sed 's/^/1\t/g' \
	| sed  -e 's/ /\t/g' -e "s/$/\t!{sample}/g" >> !{sample}.21.histo.tsv

    sed  -e 's/ /\t/g' -e "s/$/\t!{sample}/g" !{sample}.21.histo.tmp.tsv >> !{sample}.21.histo.tsv
    sed -i '/\t0\t/d' !{sample}.21.histo.tsv


    echo -e "FREQUENCY\tCOUNT\tSAMPLE" > !{sample}.13.histo.tsv
    cat !{sample}.13.kmc.json | jq -r '.Stats."#k-mers_below_min_threshold"' \
	| sed 's/^/1\t/g' \
	| sed  -e 's/ /\t/g' -e "s/$/\t!{sample}/g" >> !{sample}.13.histo.tsv

    sed  -e 's/ /\t/g' -e "s/$/\t!{sample}/g" !{sample}.13.histo.tmp.tsv >> !{sample}.13.histo.tsv
    sed -i '/\t0\t/d' !{sample}.13.histo.tsv
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

             readsPair | join(readsSingle) | (pNonpareil & pKMC)

      emit:
        nonpareilIndex = pNonpareil.out.nonpareilIndex
        kmerFrequencies = pKMC.out.histogram 
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
