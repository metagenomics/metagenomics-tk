nextflow.enable.dsl=2

include { wSaveSettingsFile } from './modules/config/module'
include { wQualityControlFile } from './modules/qualityControl/module'
include { wAssemblyFile; wAssemblyList } from './modules/assembly/module'
include { wBinning } from './modules/binning/module.nf'
include { wMagAttributesFile; wMagAttributesList; wCMSeqWorkflowFile; } from './modules/magAttributes/module.nf'
include { wDereplicateFile; wDereplicateList} from './modules/dereplication/pasolli/module'
include { wListReadMappingBwa; wFileReadMappingBwa} from './modules/readMapping/bwa/module'
include { wAnalyseMetabolites } from './modules/metabolomics/module'
include { wUnmappedReadsList; wUnmappedReadsFile } from './modules/sampleAnalysis/module'
include { wFragmentRecruitmentList; wFragmentRecruitmentFile } from './modules/fragmentRecruitment/frhit/module'

include { wAnnotateFile; wAnnotateList as wAnnotateBinsList; \
	  wAnnotateList as wAnnotateUnbinnedList; wAnnotateList as wAnnotatePlasmidList } from './modules/annotation/module'

include { wCooccurrenceList; wCooccurrenceFile } from './modules/cooccurrence/module'
include { wPlasmidsList; } from './modules/plasmids/module'


def mapJoin(channel_a, channel_b, key_a, key_b){
    channel_a \
        | map{ it -> [it[key_a], it] } \
        | cross(channel_b | map{it -> [it[key_b], it]}) \
        | map { it[0][1] + it[1][1] }
}

workflow wDereplication {
   wDereplicateFile(Channel.from(file(params?.steps?.dereplication?.pasolli?.input)))
}

workflow wAssembly {
   wAssemblyFile()
}

workflow wReadMapping {
   wFileReadMappingBwa()
}

workflow wDereplicationPath {
   wDereplicatePath()
}

workflow wCMSeqWorfklowFile {
   wCMSeqWorkflowFile(Channel.from(params?.steps?.matAttributes?.input?.genomes), Channel.from(params?.steps?.matAttributes?.input?.alignments))
}

workflow wMagAttributes {
   print(params.modules.magAttributes.getClass())
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

     wListReadMappingBwa(samplesPaired, samplesSingle, representativesList)

     binsStats | wAnalyseMetabolites

     wCooccurrenceList(wListReadMappingBwa.out.trimmedMean, gtdb)
}

workflow _wConfigurePipeline {
    
    if(params.steps.containsKey("plasmid")){
      params.steps.assembly.megahit.fastg = true
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

    wSaveSettingsFile(Channel.fromPath(params.input))

    wQualityControlFile(Channel.fromPath(params.input))

    wQualityControlFile.out.readsPair \
	| join(wQualityControlFile.out.readsSingle) | set { qcReads }

    wAssemblyList(qcReads)

    wBinning(wAssemblyList.out.contigs, qcReads)

    wBinning.out.notBinnedContigs \
	| map { notBinned -> [ notBinned[0], "notBinned", notBinned[1]]} \
	| set {notBinnedContigs}

    wBinning.out.binsStats  \
	| map{ bin -> [bin.SAMPLE, bin.BIN_ID, bin.PATH]} \
	| set { bins}

    wPlasmidsList(bins | mix(notBinnedContigs), wAssemblyList.out.fastg | join(wBinning.out.mapping))

    if(params?.steps?.fragmentRecruitment?.frhit){
       wFragmentRecruitmentList(wBinning.out.unmappedReads, Channel.fromPath(params?.steps?.fragmentRecruitment?.frhit?.genomes))
    }

    wAnnotatePlasmidList(Channel.value("meta"), wPlasmidsList.out.newPlasmids)

    wAnnotateBinsList(Channel.value("single"), bins)

    wAnnotateUnbinnedList(Channel.value("meta"), notBinnedContigs)

    wMagAttributesList(wBinning.out.bins)

    mapJoin(wMagAttributesList.out.checkm, wBinning.out.binsStats, "BIN_ID", "BIN_ID") | set { binsStats  }

    _wAggregate(wQualityControlFile.out.readsPair, wQualityControlFile.out.readsSingle, binsStats,wMagAttributesList.out.gtdb )
}
