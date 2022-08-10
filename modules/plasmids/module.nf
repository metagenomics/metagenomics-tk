nextflow.enable.dsl=2

include { pDumpLogs } from '../utils/processes'
include { pCovermContigsCoverage; pBowtie2} from '../binning/processes'

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

    tag "$sample"

    publishDir params.output, saveAs: { filename -> getOutput("${sample}", params.runid, "SCAPP", filename) }

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

    publishDir params.output, saveAs: { filename -> getOutput("${sample}", params.runid, "PLSDB", filename) }, \
            pattern: "{**.tsv}"

    containerOptions " --user 1000:1000 " +  Utils.getDockerMount(params.steps?.plasmid?.PLSDB?.database, params)

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

         _wPlasmids(samplesContigs, Channel.empty(), Channel.empty())
}


workflow wPlasmidsList {
     take:
       samplesContigs
       samplesFastg
       samplesReads
     main:
       _wPlasmids(samplesContigs, samplesFastg, samplesReads)
    emit:
      newPlasmids = _wPlasmids.out.newPlasmids
      newPlasmidsCoverage = _wPlasmids.out.newPlasmidsCoverage
}


workflow _runNonPlasmidAssemblyAnalysis {
    take:
      samplesContigs
    main:
      // Check if Contigs are plasmids
      samplesContigs | (pPlasClassLinear & pMobTyperLinear & pViralVerifyPlasmidLinear & pPlatonLinear)

      // Check which tools the user has chosen for filtering contigs
      selectedFilterTools = params?.steps?.plasmid.findAll({ tool, options -> {  options instanceof Map && options?.filter } }).collect{it.key}
      Channel.from(selectedFilterTools) | set { toolsChannel }

      // Collect output
      pPlatonLinear.out.plasmidsStats \
	| mix(pPlasClassLinear.out.probabilities) \
	| mix(pMobTyperLinear.out.plasmidsStats) \
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

      pPlasClassLinear.out.logs | mix(pMobTyperLinear.out.logs) \
	| mix(pViralVerifyPlasmidLinear.out.logs) | mix(pPlatonLinear.out.logs) | pDumpLogs

    emit:
      plasmids = samplesContigsPlasmids

}


workflow _runCircularAnalysis {
    take:
       samplesFastg
       samplesReads
    main:
       // search for new plasmids
       samplesFastg | pSCAPP
       SAMPLE_IDX = 0
       BIN_IDX = 1

       pSCAPP.out.plasmids \
	| map { plasmids -> [plasmids[SAMPLE_IDX], "plasmid_assembly", plasmids[BIN_IDX]] } \
	| set { newPlasmids }

       newPlasmids | (pPlasClassCircular & pMobTyperCircular & pViralVerifyPlasmidCircular & pPlatonCircular)

       pBowtie2(Channel.value(params.steps.plasmid?.containsKey("SCAPP")), Channel.value([Utils.getModulePath(params.modules.plasmids), \
	"SCAPP/readMapping", params.steps?.binning?.bowtie?.additionalParams?.bowtie, false]), pSCAPP.out.plasmids | join(samplesReads))

       pCovermContigsCoverage(Channel.value(true), Channel.value([Utils.getModulePath(params?.modules?.plasmids) \
	,"SCAPP/coverage", params?.steps?.plasmid?.SCAPP?.additionalParams?.coverm]), pBowtie2.out.mappedReads) 

       // Check which tools the user has chosen for filtering contigs
       selectedFilterTools = params?.steps?.plasmid.findAll({ tool, options -> {  options instanceof Map && options?.filter } }).collect{it.key}
       Channel.from(selectedFilterTools) | set { toolsChannel }

       // Collect output
       pPlatonCircular.out.plasmidsStats \
 	| mix(pPlasClassCircular.out.probabilities) \
	| mix(pMobTyperCircular.out.plasmidsStats) \
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

      pSCAPP.out.logs | mix(pPlasClassCircular.out.logs) | mix(pMobTyperCircular.out.logs)  \
	| mix(pViralVerifyPlasmidCircular.out.logs) | mix(pPlatonCircular.out.logs) | pDumpLogs

    emit:
      plasmids = filteredPlasmids
      coverage = pCovermContigsCoverage.out.coverage
}



workflow _wPlasmids {
     take:
       samplesContigs
       samplesFastg
       samplesReads
     main:
       _runNonPlasmidAssemblyAnalysis(samplesContigs)

       _runCircularAnalysis(samplesFastg, samplesReads)
    
       _runNonPlasmidAssemblyAnalysis.out.plasmids \
	| mix(_runCircularAnalysis.out.plasmids) | set { allPlasmids} 

       allPlasmids | pPLSDB

       pPLSDB.out.logs | pDumpLogs
     emit:
       newPlasmids = _runCircularAnalysis.out.plasmids
       newPlasmidsCoverage = _runCircularAnalysis.out.coverage
}
