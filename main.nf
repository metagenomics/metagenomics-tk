nextflow.enable.dsl=2
import java.util.regex.*;

include { wSaveSettingsList } from './modules/config/module'
include { pPublish as pPublishIllumina; pPublish as pPublishOnt; collectModuleFiles } from './modules/utils/processes'
include { _wGetCheckm; _wGetAssemblyFiles; _wGetIlluminaBinningFiles; } from './modules/utils/workflows'
include { wShortReadQualityControlFile; wShortReadQualityControlList} from './modules/qualityControl/shortReadQC'
include { wOntQualityControlFile; wOntQualityControlList} from './modules/qualityControl/ontQC'
include { wShortReadAssemblyFile; wShortReadAssemblyList; wTestMemSelection } from './modules/assembly/shortReadAssembler'
include { wOntAssemblyFile; wOntAssemblyList } from './modules/assembly/ontAssembler'
include { wShortReadBinningList } from './modules/binning/shortReadBinning'
include { wLongReadBinningList } from './modules/binning/ontBinning'
include { wEMGBList; _wExportPipeline } from './modules/export/emgb'
include { wMagAttributesFile; \
	wMagAttributesList as wMagAttributesList; \
	wMagAttributesList as wRecruitedGenomesAttributesList; \
	wCMSeqWorkflowFile; } from './modules/magAttributes/module.nf'
include { wDereplicateFile; wDereplicateList} from './modules/dereplication/bottomUpClustering/module'
include { wAnalyseMetabolitesList; wAnalyseMetabolitesFile } from './modules/metabolomics/module'
include { wListReadMappingBwa; wFileReadMappingBwa} from './modules/readMapping/mapping.nf'
include { wFragmentRecruitmentFile; wFragmentRecruitmentList;} from './modules/fragmentRecruitment/module'
include { wAnnotateFile; wAnnotateList as wAnnotateBinsList; \
          _wCreateProkkaInput; _wCreateProkkaGtdbInput;  \
	  wAnnotateList as wAnnotateUnbinnedList; \
	  wAnnotateList as wAnnotatePlasmidList; \
	  wAnnotateList as wAnnotateRecruitedGenomesList; } from './modules/annotation/module'
include { wCooccurrenceList; wCooccurrenceFile } from './modules/cooccurrence/module'
include { wPlasmidsList; wPlasmidsPath; } from './modules/plasmids/module'
include { wInputFile } from './modules/input/module'


def checkVersions(){
  supportedVersionsFormatted = params.supportedVersions.collect({ version -> version.YEAR + "." + version.MONTH })
  isVersionMatched = params.supportedVersions.collect({ version -> nextflow.version.matches(">=" + version.YEAR + "." + version.MONTH) \
	&& nextflow.version.matches("<" + version.YEAR + "." + (version.MONTH + 1)) }).any()

  if(!isVersionMatched){
    println "The meta-omics-toolkit was tested against the following Nextflow versions $supportedVersionsFormatted -- You are running version $nextflow.version"
    println "You can still use a different Nextflow version by using --skipVersionCheck "
    exit 1
  }
}

if(! params.skipVersionCheck){
  checkVersions()
}


def mapJoin(channel_a, channel_b, key_a, key_b){
    channel_a \
        | map{ it -> [it[key_a], it] } \
        | cross(channel_b | map{it -> [it[key_b], it]}) \
        | map { it[0][1] + it[1][1] }
}

workflow wDereplication {
   wSaveSettingsList(Channel.value("AGGREGATED"))
   wDereplicateFile(Channel.from(file(params?.steps?.dereplication?.bottomUpClustering?.input)))
}

workflow wShortReadAssembly {
   wShortReadAssemblyFile()
}

workflow wOntAssembly {
   wOntAssemblyFile()
}

workflow wOntQualityControl {
   wOntQualityControlFile()
}

workflow wShortReadQualityControl {
   wShortReadQualityControlFile()
}

workflow wReadMapping {
   wFileReadMappingBwa()
}

workflow wGetModuleVersions {
        params.modules.each { println "$it.value.name: $it.value.version.major" \
       + "." + "$it.value.version.minor" \
       + "." + "$it.value.version.patch" }
}

workflow wSRATable {
   SAMPLE_IDX = 0
   FASTQ_FILE_LEFT_IDX = 2 
   FASTQ_FILE_RIGHT_IDX = 3 
   INSTRUMENT_IDX = 1

   wInputFile() | branch {  
        ONT: it.TYPE == "OXFORD_NANOPORE"
        ILLUMINA: it.TYPE == "ILLUMINA"
   } | set { input }
     
   input.ILLUMINA | map { sample ->  [ sample.SAMPLE, sample.TYPE, sample.READS1, sample.READS2 ] } \
	| collectFile(newLine: true, seed: "SAMPLE\tINSTRUMENT\tREADS1\tREADS2"){ it -> [ "samplesILLUMINA.tsv", it[SAMPLE_IDX] \
        + "\t" + it[INSTRUMENT_IDX] \
	+ "\t" + it[FASTQ_FILE_LEFT_IDX].toString() \
	+ "\t" + it[FASTQ_FILE_RIGHT_IDX].toString()] } \
	| set {illuminaFile} 
   illuminaFile | view({ it -> it.text })
   pPublishIllumina(params.logDir, illuminaFile)

   input.ONT | map { sample ->  [ sample.SAMPLE, sample.TYPE, sample.READS ] } \
	| collectFile(newLine: true, seed: "SAMPLE\tINSTRUMENT\tREADS"){ it -> [ "samplesONT.tsv", it[SAMPLE_IDX] \
        + "\t" + it[INSTRUMENT_IDX] \
	+ "\t" + it[FASTQ_FILE_LEFT_IDX].toString() ]} \
	| set {ontFile} 
   ontFile | view({ it -> it.text })
   pPublishOnt(params.logDir, ontFile)
}

workflow wPlasmids {
   wPlasmidsPath()
}

workflow wCMSeqWorfklowFile {
   wCMSeqWorkflowFile(Channel.fromPath(params?.steps?.magAttributes?.input?.genomes), Channel.fromPath(params?.steps?.magAttributes?.input?.alignments))
}

workflow wMagAttributes {
   wMagAttributesFile(Channel.fromPath(params?.steps?.magAttributes?.input))
}

workflow wFragmentRecruitment {
   wFragmentRecruitmentFile()
}

workflow wAnnotate {
   wAnnotateFile(Channel.from(file(params?.steps?.annotation?.input)))
}

workflow wMetabolomics {
   wAnalyseMetabolitesFile()
}

workflow wCooccurrence {
   wSaveSettingsList(Channel.value("AGGREGATED"))
   wCooccurrenceFile()
}

/*
* This method either returns file path or url
*/
def getPath(f){
  return params.input.startsWith("s3://")? "s3:/" + f: f
}

workflow _wFindSamplesONT { 
   take:
     sraFiles
   main:
    // get Nanopore Fastq files
    Pattern ontFastqPattern = Pattern.compile('.*/qcONT/' + params.modules.qcONT.version.major + '..*/.*/.*qc.fq.gz$')
    sraFiles | filter({ sra, path -> ontFastqPattern.matcher(path.toString()).matches()}) \
     | map{ sra,f -> [sra, getPath(f)] } | set { ontSamples }

    // get Nanopore median phred quality score
    SAMPLE_NAME=0
    SAMPLE_STATS=1

    Pattern ontStatsPattern = Pattern.compile('.*/qcONT/' + params.modules.qcONT.version.major + '..*/.*/NanoStats.tsv$')
    sraFiles | filter({ sra, path -> ontStatsPattern.matcher(path.toString()).matches()}) \
     | map{ sra,f -> [sra, getPath(f)] } | splitCsv(sep: '\t', header: true) \
     | map { sample -> [sample[SAMPLE_NAME], sample[SAMPLE_STATS].median_qual]} |set { ontMedianQuality }

    Pattern binsONTPattern = Pattern.compile('.*/binningONT/' + params.modules.binningONT.version.major + '..*/.*/.*_bin.*.fa$')
    sraFiles | filter({ sra, path -> binsONTPattern.matcher(path.toString()).matches()}) \
     | map{ sra,f -> [SAMPLE:sra, PATH: getPath(f), BIN_ID:file(f).name] } \
     | set{ ontBins }

    // get ont binning stats
    Pattern ontBinsStatsPattern = Pattern.compile('.*/binningONT/' + params.modules.binningONT.version.major + '..*/.*/.*_bins_stats.tsv$')
    sraFiles | filter({ sra, path -> ontBinsStatsPattern.matcher(path.toString()).matches()}) \
     | splitCsv(header: true, sep: '\t') | map { sra, bins -> bins } | set{ ontBinStats }

   emit:
     ontSamples = ontSamples
     ontMedianQuality = ontMedianQuality
     ontBins = ontBins
     ontBinStats = ontBinStats
}


workflow _wFindSamplesIllumina { 
    take:
      binningFiles
      qcFiles
    main:
      binningFiles | _wGetIlluminaBinningFiles \
        | mix(qcFiles) | set { sraFiles }

      // get Illumina paired Fastq files
      Pattern illuminaPattern = Pattern.compile('.*/qc/' + params.modules.qc.version.major + '..*/.*/.*interleaved.qc.fq.gz$')
      sraFiles | filter({ sra, path -> illuminaPattern.matcher(path.toString()).matches()}) \
       | map{ sra,f -> [sra, getPath(f)] } | set { illuminaSamples }

      // get Illumina unpaired Fastq files
      Pattern unpairedIlluminaPattern = Pattern.compile('.*/qc/' + params.modules.qc.version.major + '..*/.*/.*unpaired.qc.fq.gz$')
      sraFiles | filter({ sra, path -> unpairedIlluminaPattern.matcher(path.toString()).matches()}) \
       | map{ sra,f -> [sra, getPath(f)] } | set { unpairedIlluminaSamples }

      // get Bins
      Pattern binsIlluminaPattern = Pattern.compile('.*/binning/' + params.modules.binning.version.major + '..*/.*/.*_bin.*.fa$')
      sraFiles | filter({ sra, path -> binsIlluminaPattern.matcher(path.toString()).matches()}) \
       | map{ sra,f -> [SAMPLE:sra, PATH: getPath(f), BIN_ID:file(f).name] } \
       | set{ illuminaBins }

      // get binning stats of illumina samples
      Pattern illuminaBinsStatsPattern = Pattern.compile('.*/binning/' + params.modules.binning.version.major + '..*/.*/.*_bins_stats.tsv$')
      sraFiles | filter({ sra, path -> illuminaBinsStatsPattern.matcher(path.toString()).matches()}) \
       | splitCsv(header: true, sep: '\t') | map { sra, bins -> bins } | set{illuminaBinStats}
 
    emit:
      illuminaSamples = illuminaSamples
      unpairedIlluminaSamples = unpairedIlluminaSamples
      illuminaBins = illuminaBins
      illuminaBinStats = illuminaBinStats
}


workflow _wFragmentRecruitment {
  take:
    sraFiles
  main:
    // get genomes retrieved by fragment recruitment
    Pattern recruitedGenomesPattern = Pattern.compile('.*/fragmentRecruitment/' + params.modules.fragmentRecruitment.version.major + '..*/matches/.*$')
    sraFiles | filter({ sra, path -> recruitedGenomesPattern.matcher(path.toString()).matches()}) \
     | filter({ path -> !(path ==~ /.*command.*$/)}) \
     | map{ sra,f -> [SAMPLE:sra, PATH: getPath(f), BIN_ID:file(f).name] } \
     | unique({ bin -> bin.BIN_ID})
     | set{ recruitedGenomes }

    Pattern recruitedGenomesStatsPattern = Pattern.compile('.*/fragmentRecruitment/' + params.modules.fragmentRecruitment.version.major + '..*/stats/.*_bins_stats.tsv$')
    sraFiles | filter({ sra, path -> recruitedGenomesStatsPattern.matcher(path.toString()).matches()}) \
       | splitCsv(header: true, sep: '\t') | map { sra, bins -> bins } \
       | unique({ bin -> bin.BIN_ID}) \
       | set{recruitedGenomesStats}

   emit:
     recruitedGenomes = recruitedGenomes
     recruitedGenomesStats = recruitedGenomesStats
}

workflow _wGetSamples() {
  main:
    def runID = params.runid
    def input = params.input

    // List all available SRAIDs
    Channel.from(file(input).list()) | filter({ path -> !(path ==~ /.*summary$/) && !(path ==~ /null$/) }) \
     | filter({ path -> !(path ==~ /.*AGGREGATED$/)}) \
     | set { sraDatasets }

    sraDatasets | map { sra ->  [sra, input + "/" + sra + "/" + runID + "/" ]} \
     | set {sraIDs}
  emit:
    sraIDs
    sraDatasets
}



/*
* This workflow entry point allows to aggregate information of different samples.
* It will perform analysis steps such as dereplication, read mapping and co-occurrence.
* The input files are automatically fetched as long as they adhere to the pipeline specification document (see documentation).
*/
workflow wAggregatePipeline {
    def input = params.input
    def runID = params.runid

    // Save config File
    wSaveSettingsList(Channel.value("AGGREGATED"))

    _wGetSamples()
    _wGetSamples.out.sraDatasets | set { sraDatasets }
    _wGetSamples.out.sraIDs | set { sraIDs }

    // List all files in sample directories
    sraIDs | flatMap { sraID, path -> collectModuleFiles(path, sraID, [params.modules.qc]) } | set { qcFiles }
    sraIDs | flatMap { sraID, path -> collectModuleFiles(path, sraID, [params.modules.binning]) } | set { binningFiles } 

    _wFindSamplesIllumina(binningFiles, qcFiles)

    sraIDs | flatMap { sraID, path -> collectModuleFiles(path, sraID, [params.modules.binningONT])} | set { binningONTFiles }
    sraIDs | flatMap { sraID, path -> collectModuleFiles(path, sraID, [params.modules.qcONT])} \
	| mix(binningONTFiles) | set { binningONT }

    _wFindSamplesONT(binningONT)

    sraIDs | flatMap { sraID, path -> collectModuleFiles(path, sraID, [params.modules.fragmentRecruitment])} |  _wFragmentRecruitment

    sraIDs | flatMap { sraID, path -> collectModuleFiles(path, sraID, [params.modules.magAttributes])} | set { selectedSRAMagAttributes}

    sraIDs | flatMap { sraID, path -> collectModuleFiles(path, sraID, [params.modules.annotation])} | set { selectedAnnotation}

    selectedSRAMagAttributes | _wGetCheckm 
    _wGetCheckm.out.checkm | set { checkm }

    // get gtdbtk summary files
    Pattern gtdbPattern = Pattern.compile('.*/magAttributes/' + params.modules.magAttributes.version.major + '..*/.*/.*_gtdbtk_combined.tsv$' )
    selectedSRAMagAttributes | filter({ sra, path -> gtdbPattern.matcher(path.toString()).matches()}) \
     | map { sraID, bins -> bins } \
     | splitCsv(sep: '\t', header: true) \
     | set { gtdb }

    // Get genome scale metabolic model files
    BIN_FILE_IDX = 0
    Pattern modelPattern = Pattern.compile('.*/metabolomics/' + params.modules.metabolomics.version.major + '..*/.*/.*model.xml$' )
    sraIDs | flatMap { sraID, path -> collectModuleFiles(path, sraID, [params.modules.metabolomics])} \
     | filter({ sra, path -> modelPattern.matcher(path.toString()).matches()}) \
     | map { sraID, model -> [sraID, model.name.split(".model.xml")[BIN_FILE_IDX], model]} \
     | set { models }

    recruitedGenomes = _wFragmentRecruitment.out.recruitedGenomes

    recruitedGenomesStats =  _wFragmentRecruitment.out.recruitedGenomesStats

    mapJoin(_wFindSamplesIllumina.out.illuminaBinStats | mix(_wFindSamplesONT.out.ontBinStats) \
        | mix(recruitedGenomesStats), checkm, "BIN_ID", "BIN_ID") \
        | set {checkmBinStats}

    mapJoin(checkmBinStats, _wFindSamplesIllumina.out.illuminaBins | mix(_wFindSamplesONT.out.ontBins) \
        | mix(recruitedGenomes), "BIN_ID", "BIN_ID") \
        | set {binsStatsComplete}

    _wFindSamplesIllumina.out.illuminaSamples | mix(_wFindSamplesONT.out.ontSamples) | view { sra, path -> "Files detected of SRA ID $sra" }

    _wAggregate(_wFindSamplesONT.out.ontSamples, _wFindSamplesONT.out.ontMedianQuality, _wFindSamplesIllumina.out.illuminaSamples, \
        _wFindSamplesIllumina.out.unpairedIlluminaSamples, binsStatsComplete, gtdb, models)
}


/*
* This workflow entry point allows to aggregate information of different samples.
* It will perform analysis steps such as dereplication, read mapping and co-occurrence.
* The input files are automatically fetched as long as they adhere to the pipeline specification document (see documentation).
*/
workflow wExportPipeline {
    _wExportPipeline()
}


workflow _wAggregate {
   take:
     samplesONT
     ontMedianQuality
     samplesPaired
     samplesSingle
     binsStats
     gtdb
     models
   main:
     representativeGenomesTempDir = params.tempdir + "/representativeGenomes"
     file(representativeGenomesTempDir).mkdirs()

     if (params?.steps?.dereplication?.useOnlyBinRefinement) {
        binsStats = binsStats.filter { it.PATH.toString().contains("magscot") || it.PATH.toString().contains("binningONT") }
     }

     wDereplicateList(binsStats)

     REPRESENTATIVES_PATH_IDX = 0
     wDereplicateList.out \
	| splitCsv(sep: '\t') \
        | map { it -> file(it[REPRESENTATIVES_PATH_IDX]) }\
        | set { representativesList }

     wListReadMappingBwa(samplesONT, ontMedianQuality, samplesPaired, samplesSingle, representativesList)

     // For the models we do not need the sample name
     wCooccurrenceList(wListReadMappingBwa.out.trimmedMeanMatrix, gtdb, models | map { model -> model.tail() })
}


/*
*
* This workflow configures the pipeline and sets additional parameters that are
* needed to fullfill the provided configuration.
*
*/
workflow _wConfigurePipeline {

    file(params.tempdir).mkdirs()

    // For plasmid detection we need the assembly graph of the assembler
    if(params.steps.containsKey("plasmid")){
       def fastg = [ fastg: true]
       params.steps.assembly.each { 
	 assembler, parameter -> params.steps.assembly.get(assembler).putAll(fastg)
       }
    }

    // If memory resources should be predicted by megahit then nonpareil and kmc
    // must be enabled
    if(params.steps?.assembly?.megahit?.resources?.RAM?.mode == "PREDICT"){
	if(!params.steps?.qc.containsKey("nonpareil")){
          def nonpareil = [ nonpareil: [additionalParams: " -v 10 -r 1234 "]]
          params.steps.qc.putAll(nonpareil) 
        }

	if(!params.steps?.qc.containsKey("kmc")){
          def kmc = [ kmc: [additionalParams: [ count: " -sm -cs10000 ", histo: " -cx50000 "]]]
          params.steps.qc.putAll(kmc) 
        }
    }
}


workflow wSaveSettings {
  inputSamples = wInputFile()
  wSaveSettingsList(inputSamples | map { it -> it.SAMPLE })
}


def flattenBins(binning){
  def chunkList = [];
  def SAMPLE_IDX = 0;
  def BIN_PATHS_IDX = 1;
  binning[BIN_PATHS_IDX].each {
     chunkList.add([binning[SAMPLE_IDX], it]);
  }
  return chunkList;
}



workflow _wProcessIllumina {
    take:
      reads
    main:
      wShortReadQualityControlList(reads)
      wShortReadQualityControlList.out.readsPair \
 	| join(wShortReadQualityControlList.out.readsSingle) | set { qcReads }
      wShortReadAssemblyList(qcReads, wShortReadQualityControlList.out.nonpareil, wShortReadQualityControlList.out.kmerFrequencies)
      wShortReadBinningList(wShortReadAssemblyList.out.contigs, qcReads)

    emit:
      contigs = wShortReadAssemblyList.out.contigs 
      notBinnedContigs = wShortReadBinningList.out.notBinnedContigs 
      bins = wShortReadBinningList.out.bins 
      binsStats = wShortReadBinningList.out.binsStats
      fastg = wShortReadAssemblyList.out.fastg
      mapping = wShortReadBinningList.out.mapping
      unmappedReads = wShortReadBinningList.out.unmappedReads
      contigCoverage = wShortReadBinningList.out.contigCoverage
      binContigMapping = wShortReadBinningList.out.binContigMapping
      readsPair = wShortReadQualityControlList.out.readsPair
      readsSingle = wShortReadQualityControlList.out.readsSingle
      readsPairSingle = qcReads
}

workflow _wProcessOnt {
    take:
      reads
    main:
      wOntQualityControlList(reads)
      wOntQualityControlList.out.reads | set { ontQCReads }
      wOntQualityControlList.out.medianQuality | set { medianQuality }
      wOntAssemblyList(ontQCReads | join(medianQuality))
      wLongReadBinningList(wOntAssemblyList.out.contigs, ontQCReads, wOntAssemblyList.out.graph, \
	 wOntAssemblyList.out.headerMapping, wOntAssemblyList.out.info, medianQuality)
    emit:
      contigs = wOntAssemblyList.out.contigs
      notBinnedContigs = wLongReadBinningList.out.notBinnedContigs 
      bins = wLongReadBinningList.out.bins 
      binsStats = wLongReadBinningList.out.binsStats
      mapping = wLongReadBinningList.out.mapping
      unmappedReads = wLongReadBinningList.out.unmappedReads
      contigCoverage = wLongReadBinningList.out.contigCoverage
      binContigMapping = wLongReadBinningList.out.binContigMapping
      reads = ontQCReads
      gfa = wOntAssemblyList.out.graph
      medianQuality = medianQuality
}


workflow _wMagAttributes {
  take:
    generatedBins
    recruitedGenomes
  main:
    wMagAttributesList(Channel.value("generated"), generatedBins)
    wRecruitedGenomesAttributesList(Channel.value("recruited"), recruitedGenomes)
 emit:
    checkm = wMagAttributesList.out.checkm | mix(wRecruitedGenomesAttributesList.out.checkm)
    gtdb = wMagAttributesList.out.gtdb | mix(wRecruitedGenomesAttributesList.out.gtdb)
    gtdbMissing = wMagAttributesList.out.gtdbMissing | mix(wRecruitedGenomesAttributesList.out.gtdbMissing)
    generatedCheckmFiles = wMagAttributesList.out.checkmFiles
    gtdbCombinedSummaryFiles = wMagAttributesList.out.gtdbCombinedSummaryFiles 
}


/*
* 
* Main workflow entrypoint. Takes list of files containing reads as input and produces assembly, binning, dereplication and metabolomics 
* results depending on the specification of the input.yml.
* 
* Input file with columns seperated by tabs: 
* Dataset_ID Left_Read Right_Read
*
* Left and right read could be https, s3 links or file path. 
*/
workflow wFullPipeline {

    _wConfigurePipeline()

    inputSamples = wInputFile()

    illumina = _wProcessIllumina(inputSamples | filter({ sample -> sample.TYPE == 'ILLUMINA' }) | map { it -> [ it.SAMPLE, it.READS1, it.READS2 ]} )

    ont = _wProcessOnt(inputSamples | filter({ sample -> sample.TYPE == 'OXFORD_NANOPORE' }) | map { it -> [ it.SAMPLE, it.READS ]} )

    ont.binsStats | mix(illumina.binsStats) | set { binsStats }

    SAMPLE_IDX = 0
    NOT_BINNED_PATH_IDX = 1
    ont.notBinnedContigs | mix(illumina.notBinnedContigs) 
       | map { notBinned -> [ notBinned[SAMPLE_IDX], notBinned[SAMPLE_IDX] + "_notBinned", notBinned[NOT_BINNED_PATH_IDX]]} \
       | set { notBinnedContigs }

    ont.binsStats | mix(illumina.binsStats) 
	| map{ bin -> [bin.SAMPLE, bin.BIN_ID, bin.PATH]} \
        | set { bins }

    illumina.fastg | set { fastg } 

    ont.gfa | set { gfa } 

    ont.mapping | mix(illumina.mapping) | set { mapping } 

    wFragmentRecruitmentList(illumina.unmappedReads, ont.unmappedReads, ont.medianQuality)

    ont.unmappedReads | mix(illumina.unmappedReads) | set { unmappedReads } 

    ont.contigCoverage | mix(illumina.contigCoverage) | set { contigCoverage } 

    ont.binContigMapping | mix(illumina.binContigMapping) | set { binContigMapping } 

    wSaveSettingsList(inputSamples | map { it -> it.SAMPLE })

    MAX_KMER = 0
    wPlasmidsList(bins | mix(notBinnedContigs), fastg | mix(gfa | combine(Channel.value(MAX_KMER))) | join(mapping) \
	,illumina.readsPairSingle, ont.reads, ont.medianQuality)

    ont.bins | mix(illumina.bins) |  set {generatedBinsFiles}

    _wMagAttributes(ont.bins | mix(illumina.bins), wFragmentRecruitmentList.out.foundGenomesPerSample)
    _wMagAttributes.out.checkm | set { checkm }
    _wMagAttributes.out.gtdb | set { gtdb }
    _wMagAttributes.out.gtdbMissing | set { gtdbMissing }
    _wMagAttributes.out.gtdbCombinedSummaryFiles | set { gtdbCombinedSummaryFiles }
    _wMagAttributes.out.generatedCheckmFiles | set { generatedCheckmFiles }

    mapJoin(checkm, binsStats | mix(wFragmentRecruitmentList.out.binsStats), "BIN_ID", "BIN_ID") \
	 |  set { binsStats }
     
    wAnnotatePlasmidList(Channel.value("plasmid"), Channel.value("meta"), \
    wPlasmidsList.out.newPlasmids | _wCreateProkkaInput, wPlasmidsList.out.newPlasmidsCoverage, \
    wPlasmidsList.out.newPlasmidsStats, wPlasmidsList.out.newPlasmids | map { [it[SAMPLE_IDX], 1] })

    generatedBinsFiles | map { sample, bins -> [sample, bins.size()] } | set { binsCounter }

    _wCreateProkkaGtdbInput(bins, gtdb, gtdbMissing)
    wAnnotateBinsList(Channel.value("binned"), Channel.value("single"), _wCreateProkkaGtdbInput.out.prokkaInput, \
        contigCoverage, binContigMapping, binsCounter)

    wFragmentRecruitmentList.out.foundGenomesPerSample | map { sample, genomes -> [sample, genomes.size()] } | set { recruitedGenomesCounter }
    wAnnotateRecruitedGenomesList(Channel.value("binned"), Channel.value("single"), wFragmentRecruitmentList.out.foundGenomesSeperated | _wCreateProkkaInput, \
    wFragmentRecruitmentList.out.contigCoverage, wFragmentRecruitmentList.out.genomeContigMapping, recruitedGenomesCounter)

    wAnnotateUnbinnedList(Channel.value("unbinned"), Channel.value("meta"), notBinnedContigs | _wCreateProkkaInput, \
    contigCoverage, binContigMapping, notBinnedContigs | map { it -> [it[SAMPLE_IDX], 1] })

    BIN_ID_IDX = 1
    PATH_IDX = 2
    wAnnotateBinsList.out.proteins \
	| map{ it -> [SAMPLE: it[SAMPLE_IDX], BIN_ID: it[BIN_ID_IDX], PROTEINS: it[PATH_IDX]] } | set { proteins}

    wAnalyseMetabolitesList(binsStats, mapJoin(checkm, proteins, "BIN_ID", "BIN_ID"))

    // If MMseqs was also executed on MAGs then the groupTuple operator should be adjusted
    // to wait for the correct number of MAGs. 
    ADDITIONAL_NOT_BINNED = 1
    MMSEQS_TAX_RUN_ON_MAGS = false 
    MMSEQS_TAX_RUN_ONLY_ON_NOT_BINNED = 1
    if(params.steps.containsKey("annotation") \
	&& params.steps.annotation.containsKey("mmseqs2_taxonomy") \
	&& params.steps.annotation.mmseqs2_taxonomy.containsKey("runOnMAGs")){
        MMSEQS_TAX_RUN_ON_MAGS = true 
    }


    // Export json files for EMGB
    wEMGBList(illumina.contigs | mix(ont.contigs),\
	mapping, \
        generatedBinsFiles, \
        gtdbCombinedSummaryFiles, \
        generatedCheckmFiles, \
        wAnnotateBinsList.out.gff | mix(wAnnotateUnbinnedList.out.gff) \
		| combine(binsCounter | map { sample, counter -> [sample, counter+ADDITIONAL_NOT_BINNED] }, by: SAMPLE_IDX), \
        wAnnotateBinsList.out.ffn | mix(wAnnotateUnbinnedList.out.ffn) 
		| combine(binsCounter | map { sample, counter -> [sample, counter+ADDITIONAL_NOT_BINNED] }, by: SAMPLE_IDX), \
        wAnnotateBinsList.out.faa | mix(wAnnotateUnbinnedList.out.faa) 
		| combine(binsCounter | map { sample, counter -> [sample, counter+ADDITIONAL_NOT_BINNED] }, by: SAMPLE_IDX), \
        wAnnotateBinsList.out.mmseqs2_taxonomy | mix(wAnnotateUnbinnedList.out.mmseqs2_taxonomy) \
		| combine(binsCounter \
		| map { sample, counter -> [sample, MMSEQS_TAX_RUN_ON_MAGS ? counter+ADDITIONAL_NOT_BINNED:MMSEQS_TAX_RUN_ONLY_ON_NOT_BINNED ] }, by: SAMPLE_IDX), \
        wAnnotateBinsList.out.mmseqs2_blast | mix(wAnnotateUnbinnedList.out.mmseqs2_blast)\
                | map { sample, type, db, blastResult -> [sample, db, blastResult] } \
		| combine(binsCounter | map { sample, counter -> [sample, counter+ADDITIONAL_NOT_BINNED] }, by: SAMPLE_IDX),  \
   )

   // Aggregate results of multiple samples (dereplication, read mapping, co-occurrence, ...)
   _wAggregate(ont.reads, ont.medianQuality, illumina.readsPair, illumina.readsSingle, binsStats, \
	gtdb,  wAnalyseMetabolitesList.out.models)
}
