include { wSaveSettingsList } from '../config/module'

include { pDumpLogs } from '../utils/processes'
include { pCovermContigsCoverage; pBowtie2; pMinimap2; pBwa; pBwa2} from '../binning/processes'

include { pPlaton as pPlatonCircular; \
          pPlaton as pPlatonLinear; \
          pViralVerifyPlasmid as pViralVerifyPlasmidCircular; \
          pViralVerifyPlasmid as pViralVerifyPlasmidLinear; \
          pMobTyper; \
          pPlasClass; \
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

    label 'highmemLarge'

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

    label 'highmemMedium'

    tag "$sample $binID"

    beforeScript Utils.getCreateDatabaseDirCommand("${params.polished.databases}")

    publishDir params.output, mode: "${params.publishDirMode}", \
	saveAs: { filename -> getOutput("${sample}", params.runid, "PLSDB", filename) }, \
        pattern: "{**.tsv}"

    containerOptions Utils.getDockerMount(params.steps?.plasmid?.PLSDB?.database, params)

    when params.steps.containsKey("plasmid") && params.steps.plasmid.containsKey("PLSDB")

    secret { "${S3_PLSDB_ACCESS}"!="" ? ["S3_PLSDB_ACCESS", "S3_PLSDB_SECRET"] : [] } 

    container "${params.mash_image}"

    input:
    tuple val(sample), val(binID), path(plasmids)

    output:
    tuple val("${sample}"), val("${binID}"), path("${binID}.tsv"), emit: allHits
    tuple val("${sample}"), val("${binID}"), path("${binID}_kmerThreshold_*.tsv"), emit: filteredHitsMetadata, optional: true
    tuple val("${binID}"), val("${output}"), val(params.LOG_LEVELS.INFO), file(".command.sh"), \
        file(".command.out"), file(".command.err"), file(".command.log"), emit: logs

    shell:
    output = getOutput("${sample}", params.runid, "PLSDB", "")
    S5CMD_PARAMS=params.steps?.plasmid?.PLSDB?.database?.download?.s5cmd?.params ?: ""
    S3_PLSDB_ACCESS=params?.steps?.plasmid?.PLSDB?.database?.download?.s5cmd && S5CMD_PARAMS.indexOf("--no-sign-request") == -1 ? "\$S3_PLSDB_ACCESS" : ""
    S3_PLSDB_SECRET=params?.steps?.plasmid?.PLSDB?.database?.download?.s5cmd && S5CMD_PARAMS.indexOf("--no-sign-request") == -1 ? "\$S3_PLSDB_SECRET" : ""
    template("plsdb.sh")
}


workflow wPlasmidsPath {
         Channel.from(file(params.steps.plasmid.input)) \
		| splitCsv(sep: '\t', header: true) \
		| map {it -> [ it.DATASET, it.BIN_ID, it.PATH ]} | set {samplesContigs}

         DATASET_IDX = 0
         wSaveSettingsList(samplesContigs | map { it -> it[DATASET_IDX] })
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

    tag "Sample: $sample, BinID: $binID"

    container "${params.ubuntu_image}"

    input:
    tuple val(sample), val(binID), path(plasmids)

    output:
    tuple val("${sample}"), val("${binID}"), path(plasmids), env(COUNT) 

    shell:
    '''
    COUNT=$(seqkit stats -T !{plasmids} | cut -d$'\t' -f 4 | tail -n 1)
    '''
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
    input | pCount | combine(chunkSize) | flatMap { sample -> \
	Utils.splitFilesIndex(Integer.parseInt(sample[COUNT_IDX]), sample[CHUNK_SIZE_IDX], [sample[SAMPLE_IDX], sample[BIN_ID_IDX], sample[FILE_IDX]]) } \
	| map({ sample, binID, binFile, start, end, chunkSize -> [sample, binID, binFile, start, end, chunkSize] }) | set { chunks }
  emit:
    chunks
}


def getSampleToolKey(sample){
  SAMPLE_IDX = 0
  BIN_ID_IDX = 1
  FILE_IDX = 2
  CHUNK_IDX = 3
  return ["${sample[SAMPLE_IDX]}_ttt_${sample[BIN_ID_IDX]}_ttt_${sample[FILE_IDX]}.tsv", sample[SAMPLE_IDX], sample[BIN_ID_IDX], sample[FILE_IDX], sample[CHUNK_IDX]]
}


/*
* This method tries to adjust the size of the input file for mobtyper in way that it does not run out of memory.
*/
workflow _wRunMobTyper {
   take:
     samplesContigs
   main:
      // Split input files in chunks
     _wSplit(samplesContigs, Channel.from(params.modules.plasmids.process.pMobTyper.defaults.inputSize)) | pMobTyper
     _wCollectChunks(pMobTyper.out.plasmidsStats)

   emit:
     plasmidsStats = _wCollectChunks.out.plasmidsStats
     logs = pMobTyper.out.logs

}


/*
* This method tries to adjust the size of the input file for plasClass in way that it does not run out of memory.
*/
workflow _wRunPlasClass {
   take:
     samplesContigs
   main:
      // Split input files in chunks
      _wSplit(samplesContigs, Channel.from(params.modules.plasmids.process.pPlasClass.defaults.inputSize)) | pPlasClass
      _wCollectChunks(pPlasClass.out.probabilities)

    emit:
      probabilities = _wCollectChunks.out.plasmidsStats
      logs = pPlasClass.out.logs
}


/**
*
* The pCollectFile process is designed to collect and concatenate plasmid output files.
* It uses the csvtk tool to concatenate all the .tsv files.
*
**/
process pCollectFile {

    label 'small'

    tag "Sample: $sample, Bin: $bin, Tool: $tool"

    container "${params.ubuntu_image}"

    input:
    tuple val(sample), val(bin), val(tool), path(toolOutputs)

    output:
    tuple val("${sample}"), val("${bin}"), val("${tool}"), path("*_toolOutputConcat_${tool}.tsv")

    shell:
    '''
    csvtk -t -T concat !{toolOutputs} > !{sample}_!{bin}_toolOutputConcat_!{tool}.tsv
    '''
}

/*
* This method collects the chunked outputs of tools where the inputs have been split.
*/
workflow _wCollectChunks {
  take:
    chunks
  main:
    CHUNK_PATH_IDX=3
    SAMPLE_IDX=0
    BIN_IDX=1
    TOOL_IDX=2
    chunks \
        | map { sample, type, tool, path, chunks -> \
	tuple(groupKey(sample + "_---_" + type + "_---_" + tool, chunks.toInteger()), [sample, type, tool, path]) } \
        | groupTuple() \
	| map { key, chunks -> [chunks[0][SAMPLE_IDX], chunks[0][BIN_IDX], chunks[0][TOOL_IDX], \
	chunks.stream().map{ elem -> elem[CHUNK_PATH_IDX] }.collect()] } | pCollectFile | set {statsFinal}
  emit:
    plasmidsStats = statsFinal
}


workflow _runNonPlasmidAssemblyAnalysis {
    take:
      samplesContigs
    main:
      // Check if Contigs are plasmids
      samplesContigs | (_wRunPlasClass & pViralVerifyPlasmidLinear & pPlatonLinear)

      // Check which tools the user has chosen for filtering contigs
      selectedFilterTools = params?.steps?.plasmid.findAll({ tool, options -> {  options instanceof Map && options?.filter } }).collect{it.key}
      selectedFilterToolsLength = selectedFilterTools.size() 

      Channel.from(selectedFilterTools) | set { toolsChannel }

      // Collect output
      pPlatonLinear.out.plasmidsStats \
	| mix(_wRunPlasClass.out.probabilities) \
	| mix(pViralVerifyPlasmidLinear.out.plasmidsStats) \
	| set { plasmidsStats }

      if(params?.steps?.plasmid.find{ it.key == "Filter" }?.value){
      	// Group outputs of multiple tools (e.g. Platon and Plasclass output) and use them for filtering
      	SAMPLE_ID_IDX = 0 
      	BIN_ID_IDX = 1 
      	TOOL_TYPE_IDX = 2 
      	CONTIG_HEADER_IDX = 3
      	FASTA_IDX = 4

      	plasmidsStats \
	 | filter({result ->  result[TOOL_TYPE_IDX] in selectedFilterTools }) \
	 | combine(samplesContigs, by: [SAMPLE_ID_IDX, BIN_ID_IDX]) \
         | groupTuple(by: [SAMPLE_ID_IDX, BIN_ID_IDX], size: selectedFilterToolsLength) \
         | map { result -> [result[SAMPLE_ID_IDX], result[BIN_ID_IDX], selectedFilterTools.size(), result[FASTA_IDX][0], result[CONTIG_HEADER_IDX]] } \
         | pContigsPlasmidsFilter

      	pContigsPlasmidsFilter.out.plasmids | set { samplesContigsPlasmids }
      } else {
      	samplesContigs | set { samplesContigsPlasmids }
      }

      _wRunPlasClass.out.logs \
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
	| map { plasmids -> [plasmids[SAMPLE_IDX], plasmids[SAMPLE_IDX] + "_plasmid_assembly", plasmids[BIN_IDX]] } \
	| set { newPlasmids }

       newPlasmids | (_wRunPlasClass & pViralVerifyPlasmidCircular & pPlatonCircular)

       pBowtie2(Channel.value(params.steps.containsKey("plasmid") && params.steps.plasmid?.containsKey("SCAPP") \
               && params.steps?.plasmid?.SCAPP?.additionalParams.containsKey("bowtie")), \
	Channel.value([Utils.getModulePath(params.modules.plasmids), \
	"SCAPP/readMapping/bowtie", params.steps?.plasmid?.SCAPP?.additionalParams?.bowtie, \
        params.steps?.plasmid?.SCAPP?.additionalParams?.samtoolsViewBowtie, \
	false]), pSCAPP.out.plasmids | join(illuminaReads))

       pBwa(Channel.value(params.steps.containsKey("plasmid") && params.steps.plasmid?.containsKey("SCAPP") \
		&& params.steps?.plasmid?.SCAPP?.additionalParams.containsKey("bwa")), \
	Channel.value([Utils.getModulePath(params.modules.plasmids), \
	"SCAPP/readMapping/bwa", params.steps?.plasmid?.SCAPP?.additionalParams?.bwa, \
        params.steps?.plasmid?.SCAPP?.additionalParams?.samtoolsViewBwa, \
	false]), pSCAPP.out.plasmids | join(illuminaReads))

       pBwa2(Channel.value(params.steps.containsKey("plasmid") && params.steps.plasmid?.containsKey("SCAPP") \
		&& params.steps?.plasmid?.SCAPP?.additionalParams.containsKey("bwa2")), \
	Channel.value([Utils.getModulePath(params.modules.plasmids), \
	"SCAPP/readMapping/bwa2", params.steps?.plasmid?.SCAPP?.additionalParams?.bwa2, \
        params.steps?.plasmid?.SCAPP?.additionalParams?.samtoolsViewBwa2, \
	false]), pSCAPP.out.plasmids | join(illuminaReads))

       pMinimap2(Channel.value(params.steps.containsKey("plasmid") && params?.steps?.plasmid.containsKey("SCAPP")), \
	Channel.value([Utils.getModulePath(params.modules.plasmids), \
        "SCAPP/readMapping/minimap", params.steps?.plasmid?.SCAPP?.additionalParams?.minimap, \
        params.steps?.plasmid?.SCAPP?.additionalParams?.samtoolsViewMinimap, false]), \
       pSCAPP.out.plasmids | join(ontReads))

       pBowtie2.out.mappedReads | mix(pBwa.out.mappedReads, pBwa2.out.mappedReads) \
	| combine(Channel.value(DO_NOT_ESTIMATE_IDENTITY)) \
	| mix(pMinimap2.out.mappedReads | join(ontMedianQuality, by: SAMPLE_IDX)) \
	| set { covermInput  }

       pCovermContigsCoverage(Channel.value(true), Channel.value([Utils.getModulePath(params?.modules?.plasmids) \
	,"SCAPP/coverage", params?.steps?.plasmid?.SCAPP?.additionalParams?.coverm]), covermInput) 

       // Check which tools the user has chosen for filtering contigs
       selectedFilterTools = params?.steps?.plasmid.findAll({ tool, options -> {  options instanceof Map && options?.filter } }).collect{it.key}
       selectedFilterToolsLength = selectedFilterTools.size() 

       Channel.from(selectedFilterTools) | set { toolsChannel }

       // Collect output
       pPlatonCircular.out.plasmidsStats \
 	| mix(_wRunPlasClass.out.probabilities) \
	| mix(pViralVerifyPlasmidCircular.out.plasmidsStats) \
	| set { plasmidsStats }

       if(params?.steps?.plasmid.find{ it.key == "Filter" }?.value){
      	// Group outputs of multiple tools (e.g. Platon and Plasclass output) and use them for filtering
      	 SAMPLE_ID_IDX = 0 
      	 BIN_ID_IDX = 1 
      	 TOOL_TYPE_IDX = 2 
      	 CONTIG_HEADER_IDX = 3
      	 FASTA_IDX = 4

      	 plasmidsStats \
	  | filter({result ->  result[TOOL_TYPE_IDX] in selectedFilterTools }) \
          | combine(newPlasmids, by: [SAMPLE_ID_IDX, BIN_ID_IDX] ) \
          | groupTuple(by: [SAMPLE_ID_IDX, BIN_ID_IDX], size: selectedFilterToolsLength) \
          | map { result -> [result[SAMPLE_ID_IDX], result[BIN_ID_IDX], selectedFilterTools.size(), result[FASTA_IDX][0], result[CONTIG_HEADER_IDX]] } \
          | pCircularPlasmidsFilter

      	 pCircularPlasmidsFilter.out.plasmids | set { filteredPlasmids }

      } else {
      	 newPlasmids | set { filteredPlasmids }
      }

      pSCAPP.out.logs | mix(_wRunPlasClass.out.logs)  \
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
	| mix(_runCircularAnalysis.out.plasmids) | set { allPlasmids } 

       allPlasmids | (pPLSDB & _wRunMobTyper)

       pPLSDB.out.logs | mix(_wRunMobTyper.out.logs) | pDumpLogs
     emit:
       newPlasmids = _runCircularAnalysis.out.plasmids
       newPlasmidsCoverage = _runCircularAnalysis.out.coverage
}
