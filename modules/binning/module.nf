nextflow.enable.dsl=2
params.maxbin = false
params.pBowtieMode = "mapping_reads"


MODULE="binning"
VERSION="0.1.0"
def getOutput(SAMPLE, RUNID, TOOL, filename){
    return SAMPLE + '/' + RUNID + '/' + MODULE + '/' + VERSION + '/' + TOOL + '/' + filename
}


process pGetMappingQuality {

    container 'quay.io/biocontainers/samtools:1.12--h9aed4be_1'

    tag "$sample"

    publishDir params.output, saveAs: { filename -> getOutput("${sample}", params.runid, "contigMappingBowtie", filename) }

    errorStrategy 'retry'

    label 'tiny'

    input:
    tuple val(sample), path(bam)

    output:
    tuple val("${sample}"), file("${sample}_flagstat.tsv"), emit: flagstat_raw
    tuple val("${sample}"), file("${sample}_flagstat_passed.tsv"), emit: flagstat_passed
    tuple val("${sample}"), file("${sample}_flagstat_failed.tsv"), emit: flagstat_failed

    shell:
    template 'mapping_quality.sh'
}


process pBowtie {

    container "pbelmann/bowtie2:${params.bowtie_tag}"

    label 'large'

    tag "$sample"

    publishDir params.output, saveAs: { filename -> getOutput("${sample}", params.runid, "contigMappingBowtie", filename) }

    errorStrategy 'retry'

    cpus 28

    input:
    tuple val(sample), path(contigs), path(fastqs, stageAs: 'fastq.fq.gz')

    output:
    tuple val("${sample}"), file("${sample}.bam"), optional: true, emit: mappedReads
    tuple val("${sample}"), file("${sample}_bowtie_stats.txt"), optional: true, emit: stats

    shell:
    '''
    INDEX=!{sample}.index
    bowtie2-build --threads 28 --quiet !{contigs} $INDEX 
    bowtie2 -p !{task.cpus}  --quiet --very-sensitive -x $INDEX --interleaved fastq.fq.gz 2> !{sample}_bowtie_stats.txt | samtools view -F 3584 --threads 28 -bS - | samtools sort -l 9 --threads 28 - > !{sample}.bam
    '''
}


process pMetabat {

    container "metabat/metabat:${params.metabat_tag}"

    errorStrategy 'ignore'

    tag "$sample"

    label 'large'

    publishDir params.output, saveAs: { filename -> getOutput("${sample}", params.runid, "metabat", filename) }

    when params.steps.containsKey("binning") && params.steps.binning.containsKey("metabat")

    input:
    tuple val(sample), path(contigs), path(bam)

    output:
    tuple val("${sample}"), file("${sample}_bin.*.fa"), optional: true, emit: bins
    tuple val("${sample}"), file("*.depth.txt"), optional: true, emit: metabat_depth
    tuple val("${sample}"), file("${sample}_bins_depth.tsv"), optional: true, emit: bins_depth
    tuple val("${sample}"), file("${sample}_bins_stats.tsv"), optional: true, emit: bins_stats

    shell:
    template 'metabat.sh'
}


process pMaxBin {

    container "quay.io/biocontainers/maxbin2:${params.maxbin_tag}"

  //  errorStrategy 'ignore'

    label 'large'

    tag "$sample"

    when params.maxbin

    publishDir "${params.output}/${sample}/binning/maxbin/${params.maxbin_tag}" 

    input:
    tuple val(sample), val(TYPE), path(contigs), path(reads)

    output:
    tuple val("${sample}"), env(NEW_TYPE), file("${TYPE}_${sample}/out.*.fasta"), optional: true, emit: bins

    shell:
    '''
    NEW_TYPE="!{TYPE}_maxbin"
    mkdir !{TYPE}_!{sample}
    run_MaxBin.pl -preserve_intermediate -contig !{contigs} -reads !{reads} -thread 28 -out !{TYPE}_!{sample}/out
    '''
}


def bufferBins(binning){
  def chunkList = [];
  binning[2].each {
     chunkList.add([binning[0],binning[1], it]);
  }
  return chunkList;
}


def aslist(element){
  if(element instanceof Collection){
    return element;
  } else {
    return [element];
  }
}

/*
*
* Method takes a list of the form [SAMPLE, [BIN1 path, BIN2 path]] as input
* and produces a flattend list of the form [SAMPLE, BIN 1 path, BIN 2 path]
*
*/
def flattenBins(binning){
  def chunkList = [];
  def SAMPLE_IDX = 0;
  def BIN_PATHS_IDX = 1;
  binning[BIN_PATHS_IDX].each {
     chunkList.add([binning[SAMPLE_IDX], it]);
  }
  return chunkList;
}


def mapJoin(channel_a, channel_b, key_a, key_b){
    channel_a \
        | map{ it -> [it[key_a], it] } \
        | cross(channel_b | map{it -> [it[key_b], it]}) \
        | map { it[0][1] + it[1][1] }
}

def setID(binning){
  def chunkList = [];
  binning.each {
     def bin = file(it[1]) 
     def sample = it[0]
     def binMap = [BIN_ID:bin.name, SAMPLE:sample, PATH:bin]
     chunkList.add(binMap)
  };
  return chunkList;
}

workflow wBinning {
   take: 
     contigs
     input_reads
   main:
     contigs | join(input_reads | mix(input_reads), by: 0) | pBowtie
     pBowtie.out.mappedReads | pGetMappingQuality 
     
     contigs | join(pBowtie.out.mappedReads, by: [0]) | pMetabat | set { metabat }

     metabat.bins  | map({ it -> it[1] = aslist(it[1]); it  }) | set{ bins_list }

     bins_list | map { it -> flattenBins(it) } | flatMap {it -> setID(it)} | set {binMap}

     metabat.bins_stats | map { it -> file(it[1]) } | splitCsv(sep: '\t', header: true) | set { bins_stats }

     mapJoin(bins_stats, binMap, "file", "BIN_ID") | set {binMap}

     if(params.summary){
       metabat.bins_depth | collectFile(newLine: false, keepHeader: true, storeDir: params.output + "/summary/"){ item ->
         [ "metabat_bins_depth.tsv", item[1].text  ]
       }

       metabat.bins_stats | collectFile(newLine: false, keepHeader: true, storeDir: params.output + "/summary/"){ item ->
         [ "metabat_bins_depth.tsv", item[1].text  ]
       }

       pGetMappingQuality.out.flagstat_passed | collectFile(newLine: false, keepHeader: true, storeDir: params.output + "/summary/"){ item ->
         [ "flagstat_passed.tsv", item[1].text  ]
       }

       pGetMappingQuality.out.flagstat_failed | collectFile(newLine: false, keepHeader: true, storeDir: params.output + "/summary/"){ item ->
         [ "flagstat_failed.tsv", item[1].text  ]
       }
     }

   emit:
     bins_stats = binMap
     bins = bins_list
     mapping = pBowtie.out.mappedReads
}


