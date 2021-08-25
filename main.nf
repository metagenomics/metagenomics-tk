nextflow.enable.dsl=2

//import Helper

include { wSaveSettingsFile } from './modules/config/module'
include { wQualityControlFile } from './modules/qualityControl/module'
include { wAssemblyFile; wAssemblyList } from './modules/assembly/module'
include { wBinning } from './modules/binning/module.nf'
include { wMagAttributesFile; wMagAttributesList; wCMSeqWorkflowFile; } from './modules/magAttributes/module.nf'
include { wDereplicateFile; wDereplicateList; wDereplicatePath } from './modules/dereplication/pasolli/module'
include { wReadMappingBwa } from './modules/readMapping/bwa/module'
include { wAnalyseMetabolites } from './modules/metabolomics/module'
include { wUnmappedReadsList; wUnmappedReadsFile } from './modules/sampleAnalysis/module'
include { wFragmentRecruitmentList; wFragmentRecruitmentFile } from './modules/fragmentRecruitment/frhit/module'
//include { wMashScreenFile } from './modules/fragmentRecruitment/mashScreen/module'


def mapJoin(channel_a, channel_b, key_a, key_b){
    channel_a \
        | map{ it -> [it[key_a], it] } \
        | cross(channel_b | map{it -> [it[key_b], it]}) \
        | map { it[0][1] + it[1][1] }
}

workflow wDereplication {
   wDereplicateFile(Channel.from(file(params?.steps?.dereplication?.pasolli?.input)))
}

workflow wDereplicationPath {
   wDereplicatePath()
}

workflow wCMSeqWorfklowFile {
   wCMSeqWorkflowFile(Channel.from(params?.steps?.matAttributes?.input?.genomes), Channel.from(params?.steps?.matAttributes?.input?.alignments))
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

//workflow wMashScreen {
//     wMashScreenFile(Channel.fromPath(params?.steps?.fragmentRecruitment?.mashScreen?.samples))
//}

def collectFiles(dir, sra){
   def fileList = [];
   dir.eachFileRecurse { item ->
        fileList.add([sra, item]);
  }
  return fileList;
}

workflow wAggregatePipeline {
    def baseDir = params.baseDir
    def runID = params.runid

    Channel.from(file(baseDir).list()) | filter({ path -> !(path ==~ /.*summary$/)}) \
     | filter({ path -> !(path ==~ /.*AGGREGATED$/)}) | set { sraDatasets }
    sraDatasets | map { sra ->  [sra, baseDir + "/" + sra + "/" + runID + "/" ]} \
     | set {sraIDs}

    sraIDs | flatMap { sraID -> collectFiles(file(sraID[1]), sraID[0])} | set {sraFiles}

    ///meta_out/0/2/0/test1/1/qc/0.1.0/fastp/test1_interleaved.qc.fq.gz
    sraFiles | filter({ it -> (it[1] ==~ /.*\/qc\/.*\/fastp\/.*interleaved.qc.fq.gz$/)}) \
     | map{ sra,f -> [sra, baseDir.startsWith("s3://")? "s3:/" + f: f] } | set{samples}

    samples | map { it -> "${it[0]}\t${it[1]}" } \
      | collectFile(seed: "SAMPLE\tREADS", name: 'processed_reads.txt', newLine: true) | set { samplesFile }

    sraFiles | filter({ it -> (it[1] ==~ /.*\/binning\/.*\/metabat\/.*.fa$/)}) \
     | map{ sra,f -> [SAMPLE:sra, PATH: baseDir.startsWith("s3://")? "s3:/" + f: f, BIN_ID:file(f).name] } | set{bins}

    // TODO: versionining must be defined
    sraFiles | filter({ sra, path -> (path ==~ /.*\/magAttributes\/.*\/checkm\/.*.tsv$/)}) \
     | splitCsv(header: ["PATH", "SAMPLE", "BIN_ID", "Marker lineage", "# genomes", "# markers", "# marker sets", "0", "1", "2", "3", "4", "5+", "COMPLETENESS", "CONTAMINATION", "HETEROGENEITY"], sep: '\t') \
  //   | splitCsv(header: ["SAMPLE", "BIN_ID", "Marker lineage", "# genomes", "# markers", "# marker sets", "0", "1", "2", "3", "4", "5+", "COMPLETENESS", "CONTAMINATION", "HETEROGENEITY"], sep: '\t') \
     | map { sra, bins -> bins} \
     | set { checkm }

    sraFiles | filter({ sra, path -> (path ==~ /.*\/binning\/0.1.0\/metabat\/.*_bins_stats.tsv$/)}) \
     | splitCsv(header: true, sep: '\t') | map { sra, bins -> bins } | set{binStats}

    // TODO: file vs. BIN_ID
    mapJoin(binStats, checkm, "BIN_ID", "BIN_ID") | set {checkmBinStats}
    mapJoin(checkmBinStats, bins, "file", "BIN_ID") | set {binsStatsComplete}

    _wAggregate(samplesFile, binsStatsComplete)
}

workflow _wAggregate {
   take:
     samples
     binsStats
   main:
     representativeGenomesTempDir = params.tempdir + "/representativeGenomes"
     file(representativeGenomesTempDir).mkdirs()

     wDereplicateList(binsStats)
     representativesList = wDereplicateList.out

     REPRESENTATIVES_PATH_IDX = 0
     representativesList \
       | splitCsv(sep: '\t') 
       | map { it -> file(it[REPRESENTATIVES_PATH_IDX]) } \
       | collectFile(tempDir: representativeGenomesTempDir){ item -> [ "representatives.fasta", item.text ] } \
       | set { representativesFasta }  

    // representativesList | first() | view()
     wReadMappingBwa(Channel.from('1'), representativesFasta, samples, representativesList)

     binsStats | wAnalyseMetabolites
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
   
    wSaveSettingsFile(Channel.fromPath(params.input))

    wQualityControlFile(Channel.fromPath(params.input))
    wAssemblyList(wQualityControlFile.out.processed_reads)

    wQualityControlFile.out.processed_reads \
        | map { it -> "${it[0]}\t${it[1]}" } \
        | collectFile(seed: "SAMPLE\tREADS", name: 'processed_reads.txt', newLine: true) \
        | set { samples }

    wBinning(wAssemblyList.out.contigs, wQualityControlFile.out.processed_reads)

    wUnmappedReadsList(wQualityControlFile.out.processed_reads, wBinning.out.bins)

    if(params?.steps?.fragmentRecruitment?.frhit){
       wFragmentRecruitmentList(wUnmappedReadsList.out.unmappedReads, Channel.fromPath(params?.steps?.fragmentRecruitment?.frhit?.genomes))
    }
    wMagAttributesList(wBinning.out.bins)
    mapJoin(wMagAttributesList.out.checkm, wBinning.out.bins_stats, "BIN_ID", "file") | set { binsStats  }

    _wAggregate(samples, binsStats)
}
