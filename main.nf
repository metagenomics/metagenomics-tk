nextflow.enable.dsl=2

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
include { wAnnotateS3File } from './modules/annotation/module'



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

workflow wAnnotate {
   wAnnotateS3File(Channel.from(file(params?.steps?.annotation?.input)))
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

    representativeGenomesTempDir = params.tempdir + "/representativeGenomes"
    file(representativeGenomesTempDir).mkdirs()

    wQualityControlFile(Channel.fromPath(params.input))
    wAssemblyList(wQualityControlFile.out.processed_reads)

    wQualityControlFile.out.processed_reads \
        | map { it -> "${it[0]}\t${it[1]}" } \
        | collectFile(seed: "SAMPLE\tREADS", name: 'test.txt', newLine: true) \
        | set { samples }

    wBinning(wAssemblyList.out.contigs, wQualityControlFile.out.processed_reads)

    wUnmappedReadsList(wQualityControlFile.out.processed_reads, wBinning.out.bins)

    if(params?.steps?.fragmentRecruitment?.frhit){
       wFragmentRecruitmentList(wUnmappedReadsList.out.unmappedReads, Channel.fromPath(params?.steps?.fragmentRecruitment?.frhit?.genomes))
    }
    wMagAttributesList(wBinning.out.bins)
    mapJoin(wMagAttributesList.out.checkm, wBinning.out.bins_stats, "BIN_ID", "file") | set { binsStats  }

    wDereplicateList(binsStats)
    representativesList = wDereplicateList.out

    REPRESENTATIVES_PATH_IDX = 0
    representativesList | splitCsv(sep: '\t') \
       |  map { it -> file(it[REPRESENTATIVES_PATH_IDX]) } \
       | collectFile(tempDir: representativeGenomesTempDir){ item -> [ "representatives.fasta", item.text ] } \
       | set { representativesFasta }  

    wReadMappingBwa(Channel.from('1'), representativesFasta, samples, representativesList)

    binsStats | wAnalyseMetabolites
}
