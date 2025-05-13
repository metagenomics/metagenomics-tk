nextflow.enable.dsl=2

include { pGetBinStatistics as pGetBinStatistics; \
	pCovermContigsCoverage; pCovermGenomeCoverage; pMinimap2; pMetabat } from './processes'

include { pDumpLogs } from '../utils/processes'

def getModulePath(module){
    return module.name + '/' + module.version.major + "." +
          module.version.minor + "." +
          module.version.patch
}

def getOutput(SAMPLE, RUNID, TOOL, filename){
    return SAMPLE + '/' + RUNID + '/' + getModulePath(params.modules.binningONT)  +
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


/*
*
* Binning Metagenomic Contigs via Composition, Coverage and Assembly Graphs. 
* Since MetaCoag is only used for processing Metaflye assemblies, it takes Flyes assemblyInfo as input. 
*
*/
process pMetaCoAG {

    container "${params.metacoag_image}"

    containerOptions (params.apptainer ? "" : ' --user 0:0 ' )

    tag "Sample: $sample"

    label 'highmemLarge'

    publishDir params.output, mode: "${params.publishDirMode}", saveAs: { filename -> getOutput("${sample}", params.runid, "metacoag", filename) }

    when params.steps.containsKey("binningONT") && params.steps.binningONT.containsKey("metacoag")

    input:
    tuple val(sample), path(graph), path(contigs), path(bam), path(headerMapping), path(flyeAssemblyInfo)

    output:
    tuple val("${sample}"), file("${sample}_bin.*.fa"), optional: true, emit: bins
    tuple val("${sample}"), file("${sample}_notBinned.fa"), optional: true, emit: notBinned
    tuple val("${sample}"), file("${sample}_bin_contig_mapping.tsv"), optional: true, emit: binContigMapping
    tuple file(".command.sh"), file(".command.out"), file(".command.err"), file(".command.log")

    shell:
    template "metacoag.sh"
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


workflow wLongReadBinningList {
  take:
    contigs
    inputReads
    inputGraph
    headerMapping
    assemblyInfo
    medianQuality
  main:
    _wBinning(contigs, inputReads, inputGraph, headerMapping, assemblyInfo, medianQuality)
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
     inputGraph
     headerMapping
     assemblyInfo
     medianQuality
   main:
     // Map reads against assembly and retrieve mapping quality
     SAMPLE_IDX=0

     pMinimap2(Channel.value(params?.steps?.containsKey("binningONT")), Channel.value([getModulePath(params.modules.binningONT), \
     "contigMapping", params.steps?.binningONT?.minimap?.additionalParams?.minimap, \
     params.steps?.binningONT?.minimap?.additionalParams?.samtoolsView, params.steps.containsKey("fragmentRecruitment")]), \
     contigs | join(inputReads, by: SAMPLE_IDX))

     pMinimap2.out.mappedReads | set { mappedReads }

     mappedReads | pGetMappingQuality

     pCovermContigsCoverage(Channel.value(params?.steps?.binningONT.find{ it.key == "contigsCoverage"}?.value), Channel.value([getModulePath(params.modules.binningONT), \
	"contigCoverage", params?.steps?.binningONT?.contigsCoverage?.additionalParams]), mappedReads | join(medianQuality, by: SAMPLE_IDX))

     inputGraph | join(contigs, by: SAMPLE_IDX) | join(mappedReads, by: SAMPLE_IDX) \
      | join(headerMapping, by: SAMPLE_IDX)  | join(assemblyInfo, by: SAMPLE_IDX) | pMetaCoAG 

     pMetabat(Channel.value(params?.steps?.containsKey("binningONT") && params?.steps?.binningONT.containsKey("metabat")), \
      Channel.value([getModulePath(params.modules.binningONT), \
      "metabatONT", params.steps?.binningONT?.metabat?.additionalParams]), \
      contigs | join(mappedReads, by: SAMPLE_IDX) | join(medianQuality, by: SAMPLE_IDX))

     pMetabat.out.bins | mix(pMetaCoAG.out.bins) | set { bins }


     emptyFile = file(params.tempdir + "/empty")
     emptyFile.text = ""

     ALIGNMENT_INDEX = 2
     pCovermGenomeCoverage(Channel.value(params?.steps?.binningONT.find{ it.key == "genomeCoverage"}?.value), \
        Channel.value(""),
	Channel.value([getModulePath(params.modules.binningONT), \
	"genomeCoverage", params?.steps?.binningONT?.genomeCoverage?.additionalParams]), \
	mappedReads | join(bins, by: SAMPLE_IDX) \
	| map { sample -> sample.addAll(ALIGNMENT_INDEX, emptyFile); sample } | join(medianQuality, by: SAMPLE_IDX) )

     pCovermGenomeCoverage.out.logs | pDumpLogs

     pMetabat.out.notBinned | mix(pMetaCoAG.out.notBinned) | set { notBinned }

     // Ensure that in case just one bin is produced that it still is a list
     bins | map({ it -> it[1] = aslist(it[1]); it  }) | set{ binsList }

     // Flatten metabat outputs per sample and create a map with the 
     // following entries [BIN_ID:bin.name, SAMPLE:sample, PATH:bin]
     binsList | map { it -> flattenBins(it) } | flatMap {it -> createMap(it)} | set {binMap}

     // Compute bin statistcs (e.g. N50, average coverage depth, etc. ...)
     pMetabat.out.binContigMapping | join(mappedReads, by: SAMPLE_IDX) \
	| combine(Channel.from("metabatONT")) | join(pMetabat.out.bins, by: SAMPLE_IDX) \
	| set { metabatBinStatisticsInput }
     pMetaCoAG.out.binContigMapping | join(mappedReads, by: SAMPLE_IDX) \
	| combine(Channel.from("metacoag")) | join(pMetaCoAG.out.bins, by: SAMPLE_IDX) \
	| set { metacoagBinStatisticsInput }  

     metabatBinStatisticsInput | mix(metacoagBinStatisticsInput) \
	| join(medianQuality, by: SAMPLE_IDX)| set {binStatsInput}
     pGetBinStatistics(Channel.value(getModulePath(params.modules.binningONT)), binStatsInput)

     // Add bin statistics 
     pGetBinStatistics.out.binsStats | map { it -> file(it[1]) } \
	| splitCsv(sep: '\t', header: true) | set { binsStats }
     mapJoin(binsStats, binMap, "BIN_ID", "BIN_ID") | set {binMap}
   emit:
     binsStats = binMap
     bins = binsList
     mapping = mappedReads
     notBinnedContigs = notBinned
     unmappedReads = pMinimap2.out.unmappedReads
     contigCoverage = pCovermContigsCoverage.out.coverage     
}
