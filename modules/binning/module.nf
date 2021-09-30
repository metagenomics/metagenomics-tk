nextflow.enable.dsl=2

def getOutput(SAMPLE, RUNID, TOOL, filename){
    return SAMPLE + '/' + RUNID + '/' + params.modules.binning.name + '/' + 
          params.modules.binning.version.major + "." + 
          params.modules.binning.version.minor + "." + 
          params.modules.binning.version.patch +
          '/' + TOOL + '/' + filename
}


process pGetMappingQuality {

    container "${params.samtools_image}"

    tag "$sample"

    publishDir params.output, saveAs: { filename -> getOutput("${sample}", params.runid, "contigMappingQuality", filename) }

    errorStrategy 'retry'

    label 'tiny'

    input:
    tuple val(sample), path(bam)

    output:
    tuple val("${sample}"), file("${sample}_flagstat.tsv"), emit: flagstatRaw
    tuple val("${sample}"), file("${sample}_flagstat_passed.tsv"), emit: flagstatPassed
    tuple val("${sample}"), file("${sample}_flagstat_failed.tsv"), emit: flagstatFailed
    tuple file(".command.sh"), file(".command.out"), file(".command.err"), file(".command.log")

    shell:
    template 'mapping_quality.sh'
}


process pGetBinStatistics {

    container "${params.samtools_image}"

    tag "$sample"

    publishDir params.output, saveAs: { filename -> getOutput("${sample}", params.runid, "${binner}", filename) }

    errorStrategy 'retry'

    label 'tiny'

    input:
    tuple val(sample), path(binContigMapping), path(bam), val(binner), path(bins)

    output:
    tuple val("${sample}"), file("${sample}_contigs_depth.tsv"), optional: true, emit: contigsDepth
    tuple val("${sample}"), file("${sample}_bins_stats.tsv"), optional: true, emit: binsStats
    tuple file(".command.sh"), file(".command.out"), file(".command.err"), file(".command.log")

    shell:
    template 'binStats.sh'
}


process pBowtie {

    container "${params.bowtie_image}"

    label 'large'

    tag "$sample"

    publishDir params.output, saveAs: { filename -> getOutput("${sample}", params.runid, "contigMapping", filename) }

    errorStrategy 'retry'

    input:
    tuple val(sample), path(contigs), path(fastqs, stageAs: 'fastq.fq.gz')

    output:
    tuple val("${sample}"), file("${sample}.bam"), optional: true, emit: mappedReads
    tuple val("${sample}"), file("${sample}_bowtie_stats.txt"), optional: true, emit: stats
    tuple file(".command.sh"), file(".command.out"), file(".command.err"), file(".command.log")

    shell:
    '''
    INDEX=!{sample}.index
    bowtie2-build --threads !{task.cpus} --quiet !{contigs} $INDEX 
    bowtie2 -p !{task.cpus} !{params.steps.binning.bowtie.additionalParams.bowtie} -x $INDEX --interleaved fastq.fq.gz 2> !{sample}_bowtie_stats.txt \
          | samtools view !{params.steps.binning.bowtie.additionalParams.samtoolsView} --threads !{task.cpus} -bS - \
          | samtools sort -l 9 --threads !{task.cpus} - > !{sample}.bam
    '''
}


process pMetabat {

    container "${params.metabat_image}"

    tag "$sample"

    label 'large'

    publishDir params.output, saveAs: { filename -> getOutput("${sample}", params.runid, "metabat", filename) }

    when params.steps.containsKey("binning") && params.steps.binning.containsKey("metabat")

    input:
    tuple val(sample), path(contigs), path(bam)

    output:
    tuple val("${sample}"), file("${sample}_bin.*.fa"), optional: true, emit: bins
    tuple val("${sample}"), file("${sample}_bin_contig_mapping.tsv"), optional: true, emit: binContigMapping
    tuple file(".command.sh"), file(".command.out"), file(".command.err"), file(".command.log")


    shell:
    template 'metabat.sh'
}

process pMetabinner {

    container "${params.metabinner_image}"

    tag "$sample"

    label 'large'

    publishDir params.output, saveAs: { filename -> getOutput("${sample}", params.runid, "metabinner", filename) }

    when params.steps.containsKey("binning") && params.steps.binning.containsKey("metabinner")

    containerOptions ' --user 0:0 '

    input:
    tuple val(sample), path(contigs), path(bam)

    output:
    tuple val("${sample}"), file("${sample}_bin.*.fa"), optional: true, emit: bins
    tuple val("${sample}"), file("${sample}_bin_contig_mapping.tsv"), optional: true, emit: binContigMapping
    tuple file(".command.sh"), file(".command.out"), file(".command.err"), file(".command.log")

    shell:
    template 'metabinner.sh'
}



process pMaxBin {

    container "${params.maxbin_image}"

    label 'large'

    tag "$sample"

    when params.maxbin

    publishDir "${params.output}/${sample}/binning/maxbin/${params.maxbin_tag}" 

    input:
    tuple val(sample), val(TYPE), path(contigs), path(reads)

    output:
    tuple val("${sample}"), env(NEW_TYPE), file("${TYPE}_${sample}/out.*.fasta"), optional: true, emit: bins
    tuple file(".command.sh"), file(".command.out"), file(".command.err"), file(".command.log")

    shell:
    '''
    NEW_TYPE="!{TYPE}_maxbin"
    mkdir !{TYPE}_!{sample}
    run_MaxBin.pl -preserve_intermediate -contig !{contigs} -reads !{reads} -thread 28 -out !{TYPE}_!{sample}/out
    '''
}


/*
*
* Input element is returns as an entry in a list.
*
*/
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


/*
*
* Method takes two channels with map entries and two keys as input.
* Channels are joined by the keys provided.
* Resulting channel is returned as output.
*
*/
def mapJoin(channel_a, channel_b, key_a, key_b){
    channel_a \
        | map{ it -> [it[key_a], it] } \
        | cross(channel_b | map{it -> [it[key_b], it]}) \
        | map { it[0][1] + it[1][1] }
}


/*
*
* Method takes a list of lists of the form [[SAMPLE, BIN 1 path]] 
* and produces a map of the form [BIN_ID:bin.name, SAMPLE:sample, PATH:bin]
*
*/
def createMap(binning){
  def chunkList = [];
  binning.each {
     def sample = it[0]
     def bin = file(it[1]) 
     def binMap = [BIN_ID:bin.name, SAMPLE:sample, PATH:bin]
     chunkList.add(binMap)
  };
  return chunkList;
}

workflow wBinning {
   take: 
     contigs
     inputReads
   main:
     // Map reads against assembly and retrieve mapping quality
     SAMPLE_IDX=0
     contigs | join(inputReads, by: SAMPLE_IDX) | pBowtie
     pBowtie.out.mappedReads | pGetMappingQuality 

     // Run binning tool
     contigs | join(pBowtie.out.mappedReads, by: SAMPLE_IDX) | (pMetabinner & pMetabat )
     pMetabinner.out.bins | mix(pMetabat.out.bins) | set { bins }

     // Ensure that in case just one bin is produced that it still is a list
     bins | map({ it -> it[1] = aslist(it[1]); it  }) | set{ binsList }

     // Flatten metabat outputs per sample and create a map with the 
     // following entries [BIN_ID:bin.name, SAMPLE:sample, PATH:bin]
     binsList | map { it -> flattenBins(it) } | flatMap {it -> createMap(it)} | set {binMap}

     // Compute bin statistcs (e.g. N50, average coverage depth, etc. ...)
     pMetabinner.out.binContigMapping | join(pBowtie.out.mappedReads, by: SAMPLE_IDX) \
	| combine(Channel.from("metabinner")) | join(pMetabinner.out.bins, by: SAMPLE_IDX) \
	| set { metabinnerBinStatisticsInput }  
     pMetabat.out.binContigMapping | join(pBowtie.out.mappedReads, by: SAMPLE_IDX) \
	| combine(Channel.from("metabat")) | join(pMetabat.out.bins, by: SAMPLE_IDX) \
	| set { metabatBinStatisticsInput }

     metabatBinStatisticsInput | mix(metabinnerBinStatisticsInput) | pGetBinStatistics 

     // Add bin statistics 
     pGetBinStatistics.out.binsStats | map { it -> file(it[1]) } \
	| splitCsv(sep: '\t', header: true) | set { binsStats }
     mapJoin(binsStats, binMap, "BIN_ID", "BIN_ID") | set {binMap}

     // Create summary if requested
     if(params.summary){
       pGetBinStatistics.out.binsDepth \
	| collectFile(newLine: false, keepHeader: true, storeDir: params.output + "/summary/"){ item ->
         [ "metabat_bins_depth.tsv", item[1].text  ]
       }

       pGetBinStatistics.out.binsStats \
	| collectFile(newLine: false, keepHeader: true, storeDir: params.output + "/summary/"){ item ->
         [ "metabat_bins_depth.tsv", item[1].text  ]
       }

       pGetMappingQuality.out.flagstatPassed \
	| collectFile(newLine: false, keepHeader: true, storeDir: params.output + "/summary/"){ item ->
         [ "flagstat_passed.tsv", item[1].text  ]
       }

       pGetMappingQuality.out.flagstatFailed \
	| collectFile(newLine: false, keepHeader: true, storeDir: params.output + "/summary/"){ item ->
         [ "flagstat_failed.tsv", item[1].text  ]
       }
     }
   emit:
     binsStats = binMap
     bins = binsList
     mapping = pBowtie.out.mappedReads
}
