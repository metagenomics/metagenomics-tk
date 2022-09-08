nextflow.enable.dsl=2

include { pGetBinStatistics as pGetBinStatistics; \
	pGetBinStatistics as pGetNotBinnedStatistics; \
	pCovermContigsCoverage; pBowtie2; pMetabat } from './processes'

def getModulePath(module){
    return module.name + '/' + module.version.major + "." +
          module.version.minor + "." +
          module.version.patch
}

def getOutput(SAMPLE, RUNID, TOOL, filename){
    return SAMPLE + '/' + RUNID + '/' + getModulePath(params.modules.binning)  +
          '/' + TOOL + '/' + filename
}

process pGetMappingQuality {

    container "${params.samtools_image}"

    tag "$sample"

    publishDir params.output, mode: "${params.publishDirMode}", saveAs: { filename -> getOutput("${sample}", params.runid, "readMappingQuality", filename) }

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


process pMetabinner {

    container "${params.metabinner_image}"

    tag "Sample: $sample"

    label 'large'

    publishDir params.output, mode: "${params.publishDirMode}", saveAs: { filename -> getOutput("${sample}", params.runid, "metabinner", filename) }

    when params.steps.containsKey("binning") && params.steps.binning.containsKey("metabinner")

    containerOptions ' --user 0:0 '

    input:
    tuple val(sample), path(contigs), path(bam)

    output:
    tuple val("${sample}"), file("${sample}_bin.*.fa"), optional: true, emit: bins
    tuple val("${sample}"), file("${sample}_notBinned.fa"), optional: true, emit: notBinned
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



process pGraphMB {

    container "${params.graphmb_image}"

//    containerOptions "  --user 0:0 "

    tag "Sample: $sample"

    label 'large'

    publishDir params.output, mode: "${params.publishDirMode}", saveAs: { filename -> getOutput("${sample}", params.runid, "graphmb", filename) }

    when params.steps.containsKey("binning") && params.steps.binning.containsKey("metacoag")

    input:
    tuple val(sample), path(graph), path(contigs), path(bam), path(headerMapping), path(flyeAssemblyInfo), path(reads)


    //tuple val(sample), path(graph), path(contigs), path(bam), path(headerMapping), path(flyeAssemblyInfo)
    output:
    tuple val("${sample}"), file("${sample}_bin.*.fa"), optional: true, emit: bins
    tuple val("${sample}"), file("${sample}_notBinned.fa"), optional: true, emit: notBinned
    tuple val("${sample}"), file("${sample}_bin_contig_mapping.tsv"), optional: true, emit: binContigMapping
    tuple file(".command.sh"), file(".command.out"), file(".command.err"), file(".command.log")


    shell:
    template "graphmb.sh"

    //pigz -dc !{contigs} >  contigs.fa
//    minimap2 -I 64GB -d assembly.mmi assembly2.fasta # make index
//    minimap2 -I 64GB -ax map-ont assembly.mmi <reads_file> > assembly.sam
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


workflow wShortReadBinningList {
  take:
    contigs
    inputReads     
  main:
    _wBinning(contigs, inputReads)   
  emit:
     binsStats = _wBinning.out.binsStats
     bins = _wBinning.out.bins
     mapping = _wBinning.out.mapping
     notBinnedContigs = _wBinning.out.notBinnedContigs
     unmappedReads = _wBinning.out.unmappedReads
     contigCoverage = _wBinning.out.contigCoverage
}


/*
*
* This workflow takes an input_reads channel as input with the following format [SAMPLE, READS PAIRED, READS UNPAIRED]
*
*/
workflow _wBinning {
   take: 
     contigs
     inputReads
   main:
     // Map reads against assembly and retrieve mapping quality
     SAMPLE_IDX=0
     DO_NOT_ESTIMATE_IDENTITY = "-1" 

     pBowtie2(Channel.value(params?.steps?.containsKey("binning")), Channel.value([getModulePath(params.modules.binning), \
      "contigMapping", params.steps?.binning?.bowtie?.additionalParams?.bowtie, params.steps.containsKey("fragmentRecruitment")]), \
      contigs | join(inputReads, by: SAMPLE_IDX))

     pBowtie2.out.mappedReads | set { mappedReads }

     mappedReads | pGetMappingQuality

     pCovermContigsCoverage(Channel.value(params?.steps?.binning.find{ it.key == "contigsCoverage"}?.value), Channel.value([getModulePath(params.modules.binning), \
	"contigCoverage", params?.steps?.binning?.contigsCoverage?.additionalParams]), mappedReads | join(Channel.value(DO_NOT_ESTIMATE_IDENTITY), by: SAMPLE_IDX))

     contigs | join(mappedReads, by: SAMPLE_IDX) | join(Channel.value(DO_NOT_ESTIMATE_IDENTITY), by: SAMPLE_IDX) | set { binningInput }

     pMetabinner(binningInput)

     pMetabat(Channel.value(params?.steps?.containsKey("binning") && params?.steps?.binning.containsKey("metabat")), \
      Channel.value([getModulePath(params.modules.binning), \
      "metabat", params.steps?.binning?.metabat?.additionalParams]), \
      binningInput)

     pMetabinner.out.bins | mix(pMetabat.out.bins) | set { bins }

     pMetabinner.out.notBinned | mix(pMetabat.out.notBinned) | set { notBinned }

     // Ensure that in case just one bin is produced that it still is a list
     bins | map({ it -> it[1] = aslist(it[1]); it  }) | set{ binsList }

     // Flatten metabat outputs per sample and create a map with the 
     // following entries [BIN_ID:bin.name, SAMPLE:sample, PATH:bin]
     binsList | map { it -> Utils.flattenTuple(it) } | flatMap {it -> createMap(it)} | set {binMap}

     // Compute bin statistcs (e.g. N50, average coverage depth, etc. ...)
     pMetabinner.out.binContigMapping | join(mappedReads, by: SAMPLE_IDX) \
	| combine(Channel.from("metabinner")) | join(pMetabinner.out.bins, by: SAMPLE_IDX) \
	| set { metabinnerBinStatisticsInput }  
     pMetabat.out.binContigMapping | join(mappedReads, by: SAMPLE_IDX) \
	| combine(Channel.from("metabat")) | join(pMetabat.out.bins, by: SAMPLE_IDX) \
	| set { metabatBinStatisticsInput }

     metabatBinStatisticsInput | mix(metabinnerBinStatisticsInput) \
	| join(Channel.value(DO_NOT_ESTIMATE_IDENTITY), by: SAMPLE_IDX) | set {binStatsInput}

     pGetBinStatistics(Channel.value(getModulePath(params.modules.binning)), binStatsInput)

     // Add bin statistics 
     pGetBinStatistics.out.binsStats | map { it -> file(it[1]) } \
	| splitCsv(sep: '\t', header: true) | set { binsStats }
     mapJoin(binsStats, binMap, "BIN_ID", "BIN_ID") | set {binMap}
   emit:
     binsStats = binMap
     bins = binsList
     mapping = mappedReads
     notBinnedContigs = notBinned
     unmappedReads =  pBowtie2.out.unmappedReads
     contigCoverage = pCovermContigsCoverage.out.coverage     
}
