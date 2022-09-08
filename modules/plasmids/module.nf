nextflow.enable.dsl=2

include { pDumpLogs } from '../utils/processes'
include { pCovermContigsCoverage; pBowtie2; pMinimap2} from '../binning/processes'

include { pPlaton as pPlatonCircular; \
          pPlaton as pPlatonLinear; \
          pPlasClass as pPlasClassLinear; \
          pPlasClass as pPlasClassCircular; \
          pViralVerifyPlasmid as pViralVerifyPlasmidCircular; \
          pViralVerifyPlasmid as pViralVerifyPlasmidLinear; \
          pMobTyper as pMobTyperLinear; \
          pMobTyper as pMobTyperCircular;
          pFilter as pCircularPlasmidsFilter;
          pFilter as pContigsPlasmidsFilter } from './processes'

def getBasePath(SAMPLE, RUNID, TOOL){
    return SAMPLE + '/' + RUNID + '/' + params.modules.plasmids.name + '/' +
          params.modules.plasmids.version.major + "." +
          params.modules.plasmids.version.minor + "." +
          params.modules.plasmids.version.patch +
          '/' + TOOL
}

def getOutput(SAMPLE, RUNID, TOOL, filename){
    return  getBasePath(SAMPLE, RUNID, TOOL) + '/' + filename
}

process pSCAPP {

    label 'medium'

    tag "Sample: $sample"

    publishDir params.output, mode: "${params.publishDirMode}", saveAs: { filename -> getOutput("${sample}", params.runid, "SCAPP", filename) }

    when params.steps.containsKey("plasmid") && params.steps.plasmid?.containsKey("SCAPP")

    container "${params.SCAPP_image}"

    input:
    tuple val(sample), path(assemblyGraph), val(maxKmer), path(bam)

    output:
    tuple val("${sample}"), path("${sample}_plasmids.fasta.gz"), emit: plasmids, optional: true
    tuple val("${sample}"), path("${sample}_plasmids_summary_stats.tsv"), emit: plasmidsSummaryStats, optional: true
    tuple val("${sample}"), path("${sample}_plasmids_stats.tsv"), emit: plasmidsStats, optional: true
    tuple val("${sample}"), val("${output}"), val(params.LOG_LEVELS.INFO), file(".command.sh"), \
        file(".command.out"), file(".command.err"), file(".command.log"), emit: logs


    shell:
    output = getOutput("${sample}", params.runid, "SCAPP", "")
    template("scapp.sh")
}


process pPLSDB {

    label 'medium'

    tag "$sample $binID"

    beforeScript "mkdir -p ${params.polished.databases}"

    publishDir params.output, mode: "${params.publishDirMode}", \
	saveAs: { filename -> getOutput("${sample}", params.runid, "PLSDB", filename) }, \
        pattern: "{**.tsv}"

    containerOptions Utils.getDockerMount(params.steps?.plasmid?.PLSDB?.database, params)

    when params.steps.containsKey("plasmid") && params.steps.plasmid.containsKey("PLSDB")

    container "${params.mash_image}"

    input:
    tuple val(sample), val(binID), path(plasmids)

    output:
    tuple val("${sample}"), val("${binID}"), path("${sample}_${binID}.tsv"), emit: allHits
    tuple val("${sample}"), val("${binID}"), path("${sample}_${binID}_kmerThreshold_*.tsv"), emit: filteredHitsMetadata
    tuple val("${sample}_${binID}"), val("${output}"), val(params.LOG_LEVELS.INFO), file(".command.sh"), \
        file(".command.out"), file(".command.err"), file(".command.log"), emit: logs

    shell:
    output = getOutput("${sample}", params.runid, "PLSDB", "")
    template("plsdb.sh")
}


workflow wPlasmidsPath {
         Channel.from(file(params.steps.plasmid.input)) \
		| splitCsv(sep: '\t', header: true) \
		| map {it -> [ it.DATASET, it.BIN_ID, it.PATH ]} | set {samplesContigs}

         _wPlasmids(samplesContigs, Channel.empty(), Channel.empty(), Channel.empty(), Channel.empty())
}


workflow wPlasmidsList {
     take:
       samplesContigs
       assemblyGraph
       illuminaReads
       ontReads
       ontMedianQuality
     main:
       _wPlasmids(samplesContigs, assemblyGraph, illuminaReads, ontReads, ontMedianQuality)
    emit:
      newPlasmids = _wPlasmids.out.newPlasmids
      newPlasmidsCoverage = _wPlasmids.out.newPlasmidsCoverage
}

process pCount {

    label 'tiny'

    tag "Sample: $sample"

    container "${params.ubuntu_image}"

    input:
    tuple val(sample), val(binID), path(plasmids)

    output:
    tuple val("${sample}"), val("${binID}"), path(plasmids), env(COUNT) 

    shell:
    output = getOutput("${sample}", params.runid, "SCAPP", "")
    '''
    COUNT=$(seqkit stats -T !{plasmids} | cut -d$'\t' -f 4 | tail -n 1)
    '''
}

/*
* This method takes number of entries in a input file (e.g. fata entries in multi fasta file),
* the maximum number of allowed entries per chunk and the actual input (e.g. file).
* It creates a list of indices of chunks of the input file based on the input parameters.
*/
def splitFilesIndex(seqCount, chunkSize, sample){
  def chunk=seqCount.intdiv(chunkSize)
  if(seqCount.mod(chunkSize) != 0){
      chunk = chunk + 1
  }
  def chunks = []
  for(def n : 1..chunk){
      start = (n-1) * chunkSize + 1
     
      end = n * chunkSize
    
      if(end > seqCount){
          end=seqCount
      }
      chunks.add(sample + [start, end])
  }
  return chunks
}

/*
* 
* This method splits the input file in chunks.
*
*/
workflow _wSplit {
  take:
    input
    chunkSize
  main:
    SAMPLE_IDX = 0
    BIN_ID_IDX = 1 
    FILE_IDX = 2
    COUNT_IDX=3
    CHUNK_SIZE_IDX=4
    input | pCount | combine(chunkSize) \
	| flatMap { sample -> \
	splitFilesIndex(Integer.parseInt(sample[COUNT_IDX]), sample[CHUNK_SIZE_IDX], [sample[SAMPLE_IDX], sample[BIN_ID_IDX], sample[FILE_IDX]]) } \
	|  set { chunks }
  emit:
    chunks
}


def getSampleToolKey(sample){
  SAMPLE_IDX = 0
  BIN_ID_IDX = 1
  FILE_IDX = 2
  return ["${sample[SAMPLE_IDX]}_ttt_${sample[BIN_ID_IDX]}_ttt_${sample[FILE_IDX]}", sample[SAMPLE_IDX], sample[BIN_ID_IDX], sample[FILE_IDX]]
}

/*
* This method tries to adjust the size of the input file for mobtyper in way that it does not run out of memory.
*/
workflow _wRunMobTyper {
   take:
     samplesContigs
   main:
      // Split input files in chunks
      _wSplit(samplesContigs, Channel.from(params.modules.plasmids.process.pMobTyper.defaults.inputSize)) | pMobTyperLinear

      // Create per sample and bin id a composite key
      UNIQUE_SAMPLE_KEY_IDX = 0
      pMobTyperLinear.out.plasmidsStats | map { sample -> getSampleToolKey(sample) } \
	| unique { sample -> sample[UNIQUE_SAMPLE_KEY_IDX] } | set { mobTyperStatsChunk }

      // Collect all chunks of a specific sample and bin id
      STATS_IDX = 3 
      pMobTyperLinear.out.plasmidsStats \
	| collectFile(keepHeader: true){ sample -> \
	[ getSampleToolKey(sample)[UNIQUE_SAMPLE_KEY_IDX], file(sample[STATS_IDX]).text] } \
	| map { f -> [file(f).name, f] } | set { mobTyperStatsCombined}

      // get initial ids such as sample and bin id and remove composite key
      UNIQUE_CHUNK_IDX = 0
      COMPOSED_INDEX_IDX = 0
      mobTyperStatsChunk | join(mobTyperStatsCombined, by: UNIQUE_CHUNK_IDX) \
	| map { sample -> sample.remove(COMPOSED_INDEX_IDX); sample } | set { mobTyperStatsFinal }
    emit:
      plasmidsStats = mobTyperStatsFinal
      logs = pMobTyperLinear.out.logs
}


workflow _runNonPlasmidAssemblyAnalysis {
    take:
      samplesContigs
    main:
      // Check if Contigs are plasmids
      samplesContigs | (pPlasClassLinear & _wRunMobTyper & pViralVerifyPlasmidLinear & pPlatonLinear)

      // Check which tools the user has chosen for filtering contigs
      selectedFilterTools = params?.steps?.plasmid.findAll({ tool, options -> {  options instanceof Map && options?.filter } }).collect{it.key}
      Channel.from(selectedFilterTools) | set { toolsChannel }

      // Collect output
      pPlatonLinear.out.plasmidsStats \
	| mix(pPlasClassLinear.out.probabilities) \
	| mix(_wRunMobTyper.out.plasmidsStats) \
	| mix(pViralVerifyPlasmidLinear.out.plasmidsStats) \
	| set { plasmidsStats }

      if(params?.steps?.plasmid.find{ it.key == "Filter" }?.value){
      	// Group outputs of multiple tools (e.g. Platon output and MobTyper and Plasclass) and use them for filtering
      	SAMPLE_ID_IDX = 0 
      	BIN_ID_IDX = 1 
      	TOOL_TYPE_IDX = 2 
      	CONTIG_HEADER_IDX = 3
      	FASTA_IDX = 4

      	plasmidsStats \
	 | filter({result ->  result[TOOL_TYPE_IDX] in selectedFilterTools }) \
	 | combine(samplesContigs, by: [SAMPLE_ID_IDX, BIN_ID_IDX]) \
         | groupTuple(by: [SAMPLE_ID_IDX, BIN_ID_IDX]) \
         | map { result -> [result[SAMPLE_ID_IDX], result[BIN_ID_IDX], selectedFilterTools.size(), result[FASTA_IDX][0], result[CONTIG_HEADER_IDX]] } \
         | pContigsPlasmidsFilter

      	pContigsPlasmidsFilter.out.plasmids | set { samplesContigsPlasmids }
      } else {
      	samplesContigs | set { samplesContigsPlasmids }
      }

      pPlasClassLinear.out.logs | mix(_wRunMobTyper.out.logs) \
	| mix(pViralVerifyPlasmidLinear.out.logs) | mix(pPlatonLinear.out.logs) | pDumpLogs

    emit:
      plasmids = samplesContigsPlasmids
}


workflow _runCircularAnalysis {
    take:
       assemblyGraph
       illuminaReads
       ontReads
       ontMedianQuality
    main:
       // search for new plasmids
       assemblyGraph | pSCAPP

       SAMPLE_IDX = 0
       BIN_IDX = 1
       DO_NOT_ESTIMATE_IDENTITY = "-1"

       pSCAPP.out.plasmids \
	| map { plasmids -> [plasmids[SAMPLE_IDX], "plasmid_assembly", plasmids[BIN_IDX]] } \
	| set { newPlasmids }

       newPlasmids | (pPlasClassCircular & _wRunMobTyper & pViralVerifyPlasmidCircular & pPlatonCircular)

       pBowtie2(Channel.value(params.steps.plasmid?.containsKey("SCAPP")), Channel.value([Utils.getModulePath(params.modules.plasmids), \
	"SCAPP/readMapping/bowtie", params.steps?.plasmid?.SCAPP?.additionalParams?.bowtie, false]), pSCAPP.out.plasmids | join(illuminaReads))

       pMinimap2(Channel.value(params?.steps?.plasmid.containsKey("SCAPP")), Channel.value([Utils.getModulePath(params.modules.plasmids), \
       "SCAPP/readMapping/minimap", params.steps?.plasmid?.SCAPP?.additionalParams?.minimap, false]), \
       pSCAPP.out.plasmids | join(ontReads))

       pBowtie2.out.mappedReads \
	| combine(Channel.value(DO_NOT_ESTIMATE_IDENTITY)) \
	| mix(pMinimap2.out.mappedReads | join(Channel.value(ontMedianQuality), by: SAMPLE_IDX)) \
	| set { covermInput  }

       pCovermContigsCoverage(Channel.value(true), Channel.value([Utils.getModulePath(params?.modules?.plasmids) \
	,"SCAPP/coverage", params?.steps?.plasmid?.SCAPP?.additionalParams?.coverm]), covermInput) 

       // Check which tools the user has chosen for filtering contigs
       selectedFilterTools = params?.steps?.plasmid.findAll({ tool, options -> {  options instanceof Map && options?.filter } }).collect{it.key}
       Channel.from(selectedFilterTools) | set { toolsChannel }

       // Collect output
       pPlatonCircular.out.plasmidsStats \
 	| mix(pPlasClassCircular.out.probabilities) \
	| mix(_wRunMobTyper.out.plasmidsStats) \
	| mix(pViralVerifyPlasmidCircular.out.plasmidsStats) \
	| set { plasmidsStats }

       if(params?.steps?.plasmid.find{ it.key == "Filter" }?.value){
      	// Group outputs of multiple tools (e.g. Platon output and MobTyper and Plasclass) and use them for filtering
      	 SAMPLE_ID_IDX = 0 
      	 BIN_ID_IDX = 1 
      	 TOOL_TYPE_IDX = 2 
      	 CONTIG_HEADER_IDX = 3
      	 FASTA_IDX = 4

      	 plasmidsStats \
	  | filter({result ->  result[TOOL_TYPE_IDX] in selectedFilterTools }) \
          | combine(newPlasmids, by: [SAMPLE_ID_IDX, BIN_ID_IDX] ) \
          | groupTuple(by: [SAMPLE_ID_IDX, BIN_ID_IDX]) \
          | map { result -> [result[SAMPLE_ID_IDX], result[BIN_ID_IDX], selectedFilterTools.size(), result[FASTA_IDX][0], result[CONTIG_HEADER_IDX]] } \
          | pCircularPlasmidsFilter

      	 pCircularPlasmidsFilter.out.plasmids | set { filteredPlasmids }

      } else {
      	 newPlasmids | set { filteredPlasmids }
      }

      pSCAPP.out.logs | mix(pPlasClassCircular.out.logs) | mix(_wRunMobTyper.out.logs)  \
	| mix(pViralVerifyPlasmidCircular.out.logs) | mix(pPlatonCircular.out.logs) | pDumpLogs

    emit:
      plasmids = filteredPlasmids
      coverage = pCovermContigsCoverage.out.coverage
}



workflow _wPlasmids {
     take:
       samplesContigs
       assemblyGraph
       illuminaReads
       ontReads
       ontMedianQuality
     main:
       _runNonPlasmidAssemblyAnalysis(samplesContigs)

       _runCircularAnalysis(assemblyGraph, illuminaReads, ontReads, ontMedianQuality)
    
       _runNonPlasmidAssemblyAnalysis.out.plasmids \
	| mix(_runCircularAnalysis.out.plasmids) | set { allPlasmids} 

       allPlasmids | pPLSDB

       pPLSDB.out.logs | pDumpLogs
     emit:
       newPlasmids = _runCircularAnalysis.out.plasmids
       newPlasmidsCoverage = _runCircularAnalysis.out.coverage
}
