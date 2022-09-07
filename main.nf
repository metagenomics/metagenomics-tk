nextflow.enable.dsl=2

include { wSaveSettingsList } from './modules/config/module'
include { wShortReadQualityControlFile; wShortReadQualityControlList} from './modules/qualityControl/shortReadQC'
include { wOntQualityControlFile; wOntQualityControlList} from './modules/qualityControl/ontQC'
include { wShortReadAssemblyFile; wShortReadAssemblyList } from './modules/assembly/shortReadAssembler'
include { wOntAssemblyFile; wOntAssemblyList } from './modules/assembly/ontAssembler'
include { wShortReadBinningList } from './modules/binning/shortReadBinning'
include { wLongReadBinningList } from './modules/binning/ontBinning'
include { wMagAttributesFile; wMagAttributesList; wCMSeqWorkflowFile; } from './modules/magAttributes/module.nf'
include { wDereplicateFile; wDereplicateList} from './modules/dereplication/pasolli/module'
include { wListReadMappingBwa; wFileReadMappingBwa} from './modules/readMapping/mapping.nf'
include { wAnalyseMetabolites } from './modules/metabolomics/module'
include { wUnmappedReadsList; wUnmappedReadsFile } from './modules/sampleAnalysis/module'
include { wFragmentRecruitmentList; wFragmentRecruitmentFile } from './modules/fragmentRecruitment/frhit/module'

include { wAnnotateFile; wAnnotateList as wAnnotateBinsList; \
	  wAnnotateList as wAnnotateUnbinnedList; wAnnotateList as wAnnotatePlasmidList } from './modules/annotation/module'

include { wCooccurrenceList; wCooccurrenceFile } from './modules/cooccurrence/module'
include { wPlasmidsList; wPlasmidsPath; } from './modules/plasmids/module'
include { wInputFile } from './modules/input/module'


def mapJoin(channel_a, channel_b, key_a, key_b){
    channel_a \
        | map{ it -> [it[key_a], it] } \
        | cross(channel_b | map{it -> [it[key_b], it]}) \
        | map { it[0][1] + it[1][1] }
}

workflow wDereplication {
   wDereplicateFile(Channel.from(file(params?.steps?.dereplication?.pasolli?.input)))
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
	| collectFile(newLine: true, seed: "SAMPLE\tINSTRUMENT\tREADS1\tREADS2"){ it -> [ "samplesILLUMINA", it[INSTRUMENT_IDX] \
        + "\t" + it[SAMPLE_IDX] \
	+ "\t" + it[FASTQ_FILE_LEFT_IDX].toString() \
	+ "\t" + it[FASTQ_FILE_RIGHT_IDX].toString()] } \
	| view({ it -> it.text })

   input.ONT | map { sample ->  [ sample.SAMPLE, sample.TYPE, sample.READS ] } \
	| collectFile(newLine: true, seed: "SAMPLE\tINSTRUMENT\tREADS"){ it -> [ "samplesONT", it[INSTRUMENT_IDX] \
        + "\t" + it[SAMPLE_IDX] \
	+ "\t" + it[FASTQ_FILE_LEFT_IDX].toString() ]} \
	| view({ it -> it.text })
}

workflow wDereplicationPath {
   wDereplicatePath()
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

workflow wUnmappedReads {
   wUnmappedReadsFile(Channel.fromPath(params?.steps?.sampleAnalysis?.reads), Channel.fromPath(params?.steps?.sampleAnalysis?.bins))
}

workflow wFragmentRecruitment {
   wFragmentRecruitmentFile(Channel.fromPath(params?.steps?.fragmentRecruitment?.frhit?.samples), Channel.fromPath(params?.steps?.fragmentRecruitment?.frhit?.genomes))
}

workflow wAnnotate {
   wAnnotateFile(Channel.from(file(params?.steps?.annotation?.input)))
}

workflow wCooccurrence {
   wCooccurrenceFile()
}

def collectFiles(dir, sra){
   def fileList = [];
   def moduleList = []
   params.modules.eachWithIndex { v, k -> moduleList.add(v.getKey() + "/" + v.getValue().version.major + ".") }
   
   dir.eachFileRecurse { item ->
           found = moduleList.any {  item ==~ '.*' +  it + '.*'  }
           if(found){
              fileList.add([sra, item]);
           }
   }
   return fileList;
}

workflow wAggregatePipeline {
    def baseDir = params.baseDir
    def runID = params.runid


    // List all available SRAIDs
    Channel.from(file(baseDir).list()) | filter({ path -> !(path ==~ /.*summary$/)}) \
     | filter({ path -> !(path ==~ /.*AGGREGATED$/)}) | set { sraDatasets }
    sraDatasets | map { sra ->  [sra, baseDir + "/" + sra + "/" + runID + "/" ]} \
     | set {sraIDs}

    // List all files in sample directories
    sraIDs | flatMap { sraID, path -> collectFiles(file(path), sraID)} | set {sraFiles}

    // get Fastq files
    sraFiles | filter({ sra, path -> (path ==~ /.*\/qc\/.*\/fastp\/.*interleaved.qc.fq.gz$/)}) \
     | map{ sra,f -> [sra, baseDir.startsWith("s3://")? "s3:/" + f: f] } | set{samples}

    // get available samples
    samples | map { it -> "${it[0]}\t${it[1]}" } \
      | collectFile(seed: "SAMPLE\tREADS", name: 'processed_reads.txt', newLine: true) | set { samplesFile }

    // get Bins
    sraFiles | filter({ sra, path -> (path ==~ /.*\/binning\/.*\/metabat\/.*.fa$/)}) \
     | map{ sra,f -> [SAMPLE:sra, PATH: baseDir.startsWith("s3://")? "s3:/" + f: f, BIN_ID:file(f).name] } | set{bins}

    // get Checkm results
    sraFiles | filter({ sra, path -> (path ==~ /.*\/magAttributes\/.*\/checkm\/.*.tsv$/)}) \
     | splitCsv(header: ["SAMPLE", "BIN_ID", "Marker lineage", "# genomes", "# markers", "# marker sets", "0", "1", "2", "3", "4", "5+", "COMPLETENESS", "CONTAMINATION", "HETEROGENEITY"], sep: '\t') \
     | map { sra, bins -> bins} \
     | set { checkm }

    // get gtdbtk
    sraFiles | filter({ sra, path -> (path ==~ /.*\/magAttributes\/.*\/gtdb\/.*_gtdbtk_.*.tsv$/)}) \
     | map { sra, bins -> bins} \
     | set { gtdb }

    // get binning stats
    sraFiles | filter({ sra, path -> (path ==~ /.*\/binning\/0.1.0\/metabat\/.*_bins_stats.tsv$/)}) \
     | splitCsv(header: true, sep: '\t') | map { sra, bins -> bins } | set{binStats}

    // TODO: file vs. BIN_ID
    mapJoin(binStats, checkm, "file", "BIN_ID") | set {checkmBinStats}
    mapJoin(checkmBinStats, bins, "file", "BIN_ID") | set {binsStatsComplete}

    _wAggregate(samplesFile, binsStatsComplete, gtdb)
}


workflow _wAggregate {
   take:
     samplesONT
     samplesPaired
     samplesSingle
     binsStats
     gtdb
   main:
     representativeGenomesTempDir = params.tempdir + "/representativeGenomes"
     file(representativeGenomesTempDir).mkdirs()

     wDereplicateList(binsStats)
     representativesListOfFiles = wDereplicateList.out

     REPRESENTATIVES_PATH_IDX = 0

     representativesListOfFiles \
	| splitCsv(sep: '\t') \
        | map { it -> file(it[REPRESENTATIVES_PATH_IDX]) }\
        | set { representativesList }

     wListReadMappingBwa(samplesONT, samplesPaired, samplesSingle, representativesList)

     binsStats | wAnalyseMetabolites

     wCooccurrenceList(wListReadMappingBwa.out.trimmedMean, gtdb)
}

/*
*
* This workflow configures the pipeline and sets additional parameters that are
* needed to fullfill the provided configuration.
*
*/
workflow _wConfigurePipeline {
    // For plasmid detection we need the assembly graph of the assembler
    if(params.steps.containsKey("plasmid")){
       def fastg = [ fastg: true]
       params.steps.assembly.each { 
	 assembler, parameter -> params.steps.assembly.get(assembler).putAll(fastg)
       }
    }
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
      wShortReadAssemblyList(qcReads)
      wShortReadBinningList(wShortReadAssemblyList.out.contigs, qcReads)
    emit:
      notBinnedContigs = wShortReadBinningList.out.notBinnedContigs 
      bins = wShortReadBinningList.out.bins 
      binsStats = wShortReadBinningList.out.binsStats
      fastg = wShortReadAssemblyList.out.fastg
      mapping = wShortReadBinningList.out.mapping
      unmappedReads = wShortReadBinningList.out.unmappedReads
      contigCoverage = wShortReadBinningList.out.contigCoverage
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
      notBinnedContigs = wLongReadBinningList.out.notBinnedContigs 
      bins = wLongReadBinningList.out.bins 
      binsStats = wLongReadBinningList.out.binsStats
      mapping = wLongReadBinningList.out.mapping
      unmappedReads = wLongReadBinningList.out.unmappedReads
      contigCoverage = wLongReadBinningList.out.contigCoverage
      reads = ontQCReads
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
workflow wPipeline {
   
    _wConfigurePipeline()

    inputSamples = wInputFile()

    illumina = _wProcessIllumina(inputSamples | filter({ sample -> sample.TYPE == 'ILLUMINA' }) | map { it -> [ it.SAMPLE, it.READS1, it.READS2 ]} )

    ont = _wProcessOnt(inputSamples | filter({ sample -> sample.TYPE == 'OXFORD_NANOPORE' }) | map { it -> [ it.SAMPLE, it.READS ]} )

    ont.binsStats | mix(illumina.binsStats) | set { binsStats }

    ont.notBinnedContigs | mix(illumina.notBinnedContigs) 
       | map { notBinned -> [ notBinned[0], "notBinned", notBinned[1]]} \
       | set { notBinnedContigs }

    ont.binsStats | mix(illumina.binsStats) 
	| map{ bin -> [bin.SAMPLE, bin.BIN_ID, bin.PATH]} \
        | set { bins }

    illumina.fastg | set { fastg } 

    ont.mapping | mix(illumina.mapping) | set { mapping } 

    ont.unmappedReads | mix(illumina.unmappedReads) | set { unmappedReads } 

    ont.contigCoverage | mix(illumina.contigCoverage) | set { contigCoverage } 

    wSaveSettingsList(inputSamples | map { it -> it.SAMPLE })

    wPlasmidsList(bins | mix(notBinnedContigs), fastg | join(mapping), illumina.readsPairSingle)

    if(params?.steps?.fragmentRecruitment?.frhit){
       wFragmentRecruitmentList(illumina.unmappedReads, Channel.fromPath(params?.steps?.fragmentRecruitment?.frhit?.genomes))
    }

    wMagAttributesList(ont.bins | mix(illumina.bins))

    mapJoin(wMagAttributesList.out.checkm, binsStats, "BIN_ID", "BIN_ID") | set { binsStats  }

    wAnnotatePlasmidList(Channel.value("meta"), wPlasmidsList.out.newPlasmids, null, wPlasmidsList.out.newPlasmidsCoverage)

    wAnnotateBinsList(Channel.value("single"), bins, wMagAttributesList.out.gtdb?:null, contigCoverage)

    wAnnotateUnbinnedList(Channel.value("meta"), notBinnedContigs, null, contigCoverage)

    _wAggregate(ont.reads, illumina.readsPair, illumina.readsSingle, binsStats, wMagAttributesList.out.gtdb )
}
