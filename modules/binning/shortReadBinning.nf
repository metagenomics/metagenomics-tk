nextflow.enable.dsl=2

include { pGetBinStatistics as pGetBinStatistics; \
	pCovermContigsCoverage; pCovermGenomeCoverage; pBowtie2; pMetabat; pBwa } from './processes'

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

     pBowtie2(Channel.value(params?.steps?.containsKey("binning") && params?.steps?.binning.containsKey("bowtie")), \
      Channel.value([getModulePath(params.modules.binning), \
      "contigMapping", params.steps?.binning?.bowtie?.additionalParams?.bowtie, \
      params.steps?.binning?.bowtie?.additionalParams?.samtoolsView, params.steps.containsKey("fragmentRecruitment")]), \
      contigs | join(inputReads, by: SAMPLE_IDX))

     pBwa(Channel.value(params?.steps?.containsKey("binning") && params?.steps?.binning.containsKey("bwa")), \
      Channel.value([getModulePath(params.modules.binning), \
      "contigMapping", params.steps?.binning?.bwa?.additionalParams?.bwa, \
      params.steps?.binning?.bwa?.additionalParams?.samtoolsView,
      params.steps.containsKey("fragmentRecruitment")]), \
      contigs | join(inputReads, by: SAMPLE_IDX))

     pBowtie2.out.mappedReads | mix(pBwa.out.mappedReads) | set { mappedReads }
     pBowtie2.out.unmappedReads | mix(pBwa.out.unmappedReads) | set { unmappedReads }

     mappedReads | pGetMappingQuality

     pCovermContigsCoverage(Channel.value(params?.steps?.binning.find{ it.key == "contigsCoverage"}?.value), \
	Channel.value([getModulePath(params.modules.binning), \
	"contigCoverage", params?.steps?.binning?.contigsCoverage?.additionalParams]), \
	mappedReads | combine(Channel.value(DO_NOT_ESTIMATE_IDENTITY)))

     
      contigs | join(mappedReads, by: SAMPLE_IDX) | set { binningInput }

     pMetabinner(binningInput)

     pMetabat(Channel.value(params?.steps?.containsKey("binning") && params?.steps?.binning.containsKey("metabat")), \
      Channel.value([getModulePath(params.modules.binning), \
      "metabat", params.steps?.binning?.metabat?.additionalParams]), \
      binningInput | combine(Channel.value(DO_NOT_ESTIMATE_IDENTITY)))

     pMetabinner.out.bins | mix(pMetabat.out.bins) | set { bins }

     emptyFile = file(params.tempdir + "/empty")
     emptyFile.text = ""

     ALIGNMENT_INDEX = 2
     pCovermGenomeCoverage(Channel.value(params?.steps?.binning.find{ it.key == "genomeCoverage"}?.value), \
	Channel.value([getModulePath(params.modules.binning), \
	"genomeCoverage", params?.steps?.binning?.genomeCoverage?.additionalParams]), \
	mappedReads | join(bins, by: SAMPLE_IDX) \
	| map { sample -> sample.addAll(ALIGNMENT_INDEX, emptyFile); sample } | combine(Channel.value(DO_NOT_ESTIMATE_IDENTITY)))

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
	| combine(Channel.value(DO_NOT_ESTIMATE_IDENTITY)) | set {binStatsInput}

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
     unmappedReads =  unmappedReads
     contigCoverage = pCovermContigsCoverage.out.coverage     
}
