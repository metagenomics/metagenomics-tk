nextflow.enable.dsl=2

//TODO: remove unneeded imports
include { pGetBinStatistics as pGetBinStatistics; \
	pCovermContigsCoverage; pCovermGenomeCoverage; pMinimap2; pBwa2 } from './processes'

def getOutput(SAMPLE, RUNID, TOOL, filename){
    return SAMPLE + '/' + RUNID + '/' + params.modules.binningHybrid.name + '/' +
          params.modules.binningHybrid.version.major + "." +
          params.modules.binningHybrid.version.minor + "." +
          params.modules.binningHybrid.version.patch +
          '/' + TOOL + '/' + filename
}

def getModulePath(module){
    return module.name + '/' + module.version.major + "." +
          module.version.minor + "." +
          module.version.patch
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

//combine ont and illumina reads into one file:
process pSamtoolsMerge {

    container "${params.samtools_image}"

    label 'small'

    publishDir params.output, mode: "${params.publishDirMode}", saveAs: { filename -> getOutput("${sample}", params.runid, "samtoolsMerge", filename) }

    input:
      tuple val(sample), path(ontBam, stageAs: 'ontMapping.bam'), path(illuminaBam, stageAs: 'illuminaMapping.bam')

    output:
      tuple val("${sample}"), file("${sample}_mapping.bam"), emit: combinedBAM
      tuple file(".command.sh"), file(".command.out"), file(".command.err"), file(".command.log")

    shell:
    '''
    samtools merge !{sample}_mapping.bam --threads !{task.cpus} !{ontBam} !{illuminaBam}
    '''
}

process pMetabatHybrid {

    container "${params.metabat_image}"

    label 'large'

    when params.steps.containsKey("binningHybrid") && params?.steps?.binningHybrid.containsKey("metabat")

    publishDir params.output, mode: "${params.publishDirMode}", saveAs: { filename -> getOutput("${sample}", params.runid, "metabat", filename) }

    input:
      tuple val(sample), path(contigs, stageAs: 'contigs.fasta'), path(ontBam, stageAs: 'ontMapping.bam'), path(illuminaBam, stageAs: 'illuminaMapping.bam')

    output:
      tuple val("${sample}"), file("${sample}_bin.*.fa"), optional: true, emit: bins
      tuple val("${sample}"), file("${sample}_notBinned.fa"), optional: true, emit: notBinned
      tuple val("${sample}"), file("${sample}_bin_contig_mapping.tsv"), optional: true, emit: binContigMapping
      tuple file(".command.sh"), file(".command.out"), file(".command.err"), file(".command.log")


    //TODO: estimate quality like in ONT binning?
    shell:
      template 'metabatHybrid.sh'
}

workflow wHybridBinningList {
  take:
    contigs
    combinedReads
  main:
    _wHybridBinning(contigs, combinedReads)
  emit:
    binsStats = _wHybridBinning.out.binsStats
    bins = _wHybridBinning.out.bins
    mapping = _wHybridBinning.out.mapping
    notBinnedContigs = _wHybridBinning.out.notBinnedContigs
    unmappedReads = _wHybridBinning.out.unmappedReads
    contigCoverage = _wHybridBinning.out.contigCoverage
}

workflow _wHybridBinning {
  take:
    contigs
    combinedReads
  main:
    SAMPLE_IDX = 0
    ONT_IDX = 1
    PAIRED_END_IDX = 2
    UNPAIRED_END_IDX = 3
    ONT_MEDIAN_QUALITY_IDX = 4
    
    illuminaReads = combinedReads | map {it -> [it[SAMPLE_IDX], it[PAIRED_END_IDX], it[UNPAIRED_END_IDX]] }
    ontReads = combinedReads | map {it -> [it[SAMPLE_IDX], it[ONT_IDX]] }
    medianQuality = combinedReads | map {it -> [it[SAMPLE_IDX], it[ONT_MEDIAN_QUALITY_IDX]] }
    
    //1) Run BWA2 mapping    
    bwaInput = contigs | join(illuminaReads, by: SAMPLE_IDX) 
    pBwa2(Channel.value(params?.steps?.containsKey("binningHybrid") && params?.steps?.binningHybrid.containsKey("bwa2")), \
      Channel.value([getModulePath(params.modules.binningHybrid), \
      "illuminaMapping", params.steps?.binningHybrid?.bwa2?.additionalParams?.bwa2, \
      params.steps?.binningHybrid?.bwa2?.additionalParams?.samtoolsView,
      params.steps.containsKey("fragmentRecruitment")]), \
      bwaInput )    
    pBwa2.out.mappedReads | set { mappedIlluminaReads }


    //2) Run Minimap2 mapping
    minimapInput = contigs | join(ontReads, by: SAMPLE_IDX) 
    pMinimap2(Channel.value(params?.steps?.containsKey("binningHybrid") && params?.steps?.binningHybrid.containsKey("minimap")), \
      Channel.value([getModulePath(params.modules.binningHybrid), \
      "ontMapping", params.steps?.binningHybrid?.minimap?.additionalParams?.minimap, \
      params.steps?.binningHybrid?.minimap?.additionalParams?.samtoolsView, params.steps.containsKey("fragmentRecruitment")]), \
      minimapInput )
    pMinimap2.out.mappedReads | set { mappedONTReads }

    mappedReads = mappedONTReads | join(mappedIlluminaReads)
    unmappedReads = pMinimap2.out.unmappedReads | join(pBwa2.out.unmappedReads)

    //3) combine BAMs:
    pSamtoolsMerge(mappedReads)
    combinedBAM = pSamtoolsMerge.out.combinedBAM


    //4) Run metabat
    //TODO: run with combinedBAM?
    metabatInput = contigs | join(pMinimap2.out.mappedReads) | join(pBwa2.out.mappedReads) 
    pMetabatHybrid(metabatInput)

    pMetabatHybrid.out.bins | set { bins }
    // Ensure that in case just one bin is produced that it still is a list
    bins | map({ it -> it[1] = aslist(it[1]); it  }) | set{ binsList }

    // Flatten metabat outputs per sample and create a map with the 
    // following entries [BIN_ID:bin.name, SAMPLE:sample, PATH:bin]
    binsList | map { it -> flattenBins(it) } | flatMap {it -> createMap(it)} | set {binMap}

    pMetabatHybrid.out.notBinned | set { notBinned }
 
    //5)Coverm
    covermMappingInput = mappedONTReads | join(mappedIlluminaReads) | join (medianQuality)
    pCovermContigsCoverage(Channel.value(params?.steps?.binningHybrid.find{ it.key == "contigsCoverage"}?.value), Channel.value([getModulePath(params.modules.binningHybrid), \
	"contigCoverage", params?.steps?.binningHybrid?.contigsCoverage?.additionalParams]), combinedBAM | join (medianQuality))
    
    //6) Compute bin statistics:
     pMetabatHybrid.out.binContigMapping | join(combinedBAM, by: SAMPLE_IDX) \
	| combine(Channel.from("metabatHybrid")) | join(bins, by: SAMPLE_IDX) \
	| set { metabatBinStatisticsInput }  
     metabatBinStatisticsInput | join(medianQuality, by: SAMPLE_IDX)| set {binStatsInput}
     
     pGetBinStatistics(Channel.value(getModulePath(params.modules.binningHybrid)), binStatsInput)

     // Add bin statistics 
     pGetBinStatistics.out.binsStats | map { it -> file(it[1]) } \
	| splitCsv(sep: '\t', header: true) | set { binsStats }
     mapJoin(binsStats, binMap, "BIN_ID", "BIN_ID") | set {binMap}      

  emit:
    binsStats = binMap
    bins = binsList
    mapping = combinedBAM
    notBinnedContigs = notBinned
    unmappedReads = unmappedReads
    contigCoverage = pCovermContigsCoverage.out.coverage   
}

