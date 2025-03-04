include { wSaveSettingsList } from '../config/module'

/*
* This process runs Adapter trimming (Porechop) and quality control (filtlong) in one step in order to reduce disk space consumption of the work directory.
*/
process pPorechop {

    label 'small'

    tag "$sample"

    cache 'deep'

    publishDir params.output, mode: "${params.publishDirMode}", saveAs: { filename -> Output.getOutput("${sample}", params.runid, "porechop", params.modules.qcONT, filename) }, \
	pattern: "{.command.sh,.command.out,.command.err,.command.log}"

    when params?.steps?.containsKey("qcONT") && params?.steps?.qcONT.containsKey("porechop")

    container "${params.porechop_image}"

    input:
    tuple val(sample), path(read1, stageAs: "reads.fq.gz"), val(start), val(stop), val(chunkSize)

    output:
    tuple val("${sample}"), path("${sample}_qc.fq.gz"), val("${chunkSize}"), emit: reads
    tuple file(".command.sh"), file(".command.out"), file(".command.err"), file(".command.log")

    shell:
    '''
    seqkit range -r !{start}:!{stop} !{read1} > reads_extracted.fq
    porechop -i reads_extracted.fq -o reads.porechoped.fq.gz --threads !{task.cpus} !{params.steps.qcONT.porechop.additionalParams.porechop}
    filtlong !{params.steps.qcONT.porechop.additionalParams.filtlong} reads.porechoped.fq.gz | pigz --processes !{task.cpus} > !{sample}_qc.fq.gz
    '''
}

/**
*
* This process counts entries in a fastq file. 
*
**/
process pCount {

    label 'tiny'

    tag "Sample: $sample"

    container "${params.ubuntu_image}"

    when params?.steps?.containsKey("qcONT") && params?.steps?.qcONT.containsKey("porechop")

    input:
    tuple val(sample), path(fastq)

    output:
    tuple val("${sample}"), path(fastq), env(COUNT) 

    shell:
    '''
    set -o pipefail
    COUNT=$(seqkit stats -T <(cat !{fastq}) | cut -d$'\t' -f 4 | tail -n 1)
    '''
}

/**
*
* See pCount process
*
**/
process pCountDownload {

    label 'tiny'

    tag "Sample: $sample"

    container "${params.ubuntu_image}"

    when params?.steps?.containsKey("qcONT") && params?.steps?.qcONT.containsKey("porechop")

    input:
    tuple val(sample), env(readUrl)

    output:
    tuple val("${sample}"), path(fastq), env(COUNT) 

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

    COUNT=$(seqkit stats -T <(cat reads.fq.gz) | cut -d$'\t' -f 4 | tail -n 1)
    '''
}

/**
*
* The pCollectFile process is designed to collect and concatenate fastq files.
*
**/
process pCollectFile {

    label 'tiny'

    tag "Sample: $sample"

    container "${params.ubuntu_image}"
    publishDir params.output, mode: "${params.publishDirMode}", \
	saveAs: { filename -> Output.getOutput("${sample}", params.runid, "porechop", params.modules.qcONT, filename) }, \
        pattern: "{**_qc.fq.gz}"

    input:
    tuple val(sample), path(reads, stageAs: "reads_*")

    output:
    tuple val("${sample}"), path("${sample}_qc.fq.gz"), emit: reads

    shell:
    '''
    cat !{reads} > !{sample}_qc.fq.gz
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

    publishDir params.output, mode: "${params.publishDirMode}", saveAs: { filename -> Output.getOutput("${sample}", params.runid, "filterHumanONT", params.modules.qcONT, filename) }

    when params.steps.containsKey("qcONT") && params?.steps?.qcONT.containsKey("filterHumanONT")

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
    output = Output.getOutput("${sample}", params.runid, "filterHumanONT", params.modules.qcONT, "")
    ADDITIONAL_PARAMS=params.steps?.qcONT?.filterHumanONT?.additionalParams ?: ""
    S3_filter_ACCESS=params?.steps?.qcONT?.filterHumanONT?.database?.download?.s5cmd && S5CMD_PARAMS.indexOf("--no-sign-request") == -1 ? "\$S3_filter_ACCESS" : ""
    S3_filter_SECRET=params?.steps?.qcONT?.filterHumanONT?.database?.download?.s5cmd && S5CMD_PARAMS.indexOf("--no-sign-request") == -1 ? "\$S3_filter_SECRET" : ""
    template 'filterHumanONT.sh'
}


process pNanoPlot {

    label 'small'

    tag "Sample: $sample"

    publishDir params.output, mode: "${params.publishDirMode}", saveAs: { filename -> Output.getOutput("${sample}", params.runid, "nanoplot", params.modules.qcONT, filename) }

    when params?.steps?.containsKey("qcONT") && params?.steps?.qcONT.containsKey("nanoplot")

    container "${params.nanoplot_image}"

    input:
    tuple val(sample), path(read1, stageAs: "reads.fq.gz")

    output:
    tuple val("${sample}"), path("*.html"), emit: html
    tuple val("${sample}"), path("*.png"), emit: png
    tuple val("${sample}"), path("NanoStats.tsv"), emit: tsv
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

/*
*
* This workflow splits input fastq files into chunks.
*
*/
workflow _wSplit {
  take:
    sample
    chunkSize
  main:

    SAMPLE_IDX = 0
    FASTQ_FILE = 1
    FASTQ_LENGTH_IDX = 2
    CHUNK_SIZE_IDX = 3

    sample | combine(chunkSize) \
	| flatMap { sample -> Utils.splitFilesIndex(Integer.parseInt(sample[FASTQ_LENGTH_IDX]), sample[CHUNK_SIZE_IDX], [sample[SAMPLE_IDX], sample[FASTQ_FILE]]) } \
        | set { chunks }
  emit:
    chunks = chunks
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

             // Input files can either be directly provided as input or downloaded.
             // Sequences are then split in chunks to avoid any RAM issues
             samples.noDownload | pCount | set { notDownloadedSamples }
             samples.download | pCountDownload | set { downloadedSamples }

             PORECHOP_CHUNK_SIZE_DEFAULT=450000
             _wSplit(notDownloadedSamples \
		| mix(downloadedSamples), Channel.value(params.steps?.qcONT?.porechop?.additionalParams?.chunkSize?:PORECHOP_CHUNK_SIZE_DEFAULT)) \
		| pPorechop

             // Porechop processes chunks and the resulting fastq files are collected
             pPorechop.out.reads | map { sample, fastq, chunkSize  ->  tuple(groupKey(sample, chunkSize.toInteger()), [sample, fastq])  } \
              | groupTuple() | map { sample, dataset ->  [sample ,dataset.stream().map{ elem -> elem[FASTP_FILE_IDX] }.collect()] } \
	      | pCollectFile

             pCollectFile.out.reads | set { reads }

             filteredSeqs = Channel.empty()
             if(params.steps.containsKey("qcONT") && params.steps.qcONT.containsKey("filterHumanONT")){
             	reads | pFilterHumanONT
                pFilterHumanONT.out.filteredSeqs | set{ filteredSeqs }
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
