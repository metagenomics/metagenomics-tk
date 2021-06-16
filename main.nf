nextflow.enable.dsl=2


include { wAssemblyFile } from './modules/assembly/module'
include { wBinning } from './modules/binning/module.nf'
include { wMagAttributesFile; wMagAttributesList } from './modules/magAttributes/module.nf'
include { wDereplicateFile; wDereplicateList } from './modules/dereplication/pasolli/module'
include { wReadMappingBwa } from './modules/readMapping/bwa/module'
include { wAnalyseMetabolites } from './modules/metabolomics/module'
include { wUnmappedReadsList; wUnmappedReadsFile } from './modules/sampleAnalysis/module'
include { wFragmentRecruitmentList; wFragmentRecruitmentFile } from './modules/fragmentRecruitment/frhit/module'


def mapJoin(channel_a, channel_b, key){
    channel_a \
        | map{ it -> [it[key], it] } \
        | cross(channel_b | map{it -> [it[key], it]}) \
        | map { it[0][1] + it[1][1] }
}

workflow wDereplication {
   wDereplicateFile(Channel.from(file(params?.steps?.dereplication?.pasolli?.input)))
}


workflow wMagAttributes {
   wMagAttributesFile(Channel.fromPath(params?.steps?.magAttributes?.input))
}

workflow run_bwa {
//    bwa(Channel.from('1'), Channel.from(params.input), Channel.fromPath(params.mapping_samples),Channel.fromPath(params.list_of_representatives))
}


workflow wUnmappedReads {
     wUnmappedReadsFile(Channel.fromPath(params?.steps?.sampleAnalysis?.reads), Channel.fromPath(params?.steps?.sampleAnalysis?.bins))
}


workflow wFragmentRecruitment {
     wFragmentRecruitmentFile(Channel.fromPath(params?.steps?.fragmentRecruitment?.frhit?.samples), Channel.fromPath(params?.steps?.fragmentRecruitment?.frhit?.genomes))
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
    representativeGenomesTempDir = params.tempdir + "/representativeGenomes"
    file(representativeGenomesTempDir).mkdirs()

    wAssemblyFile(Channel.fromPath(params.input))

    wAssemblyFile.out.processed_reads \
        | map { it -> "${it[0]}\t${it[1]}" } \
        | collectFile(seed: "SAMPLE\tREADS", name: 'test.txt', newLine: true) \
        | set { samples }

    wBinning(wAssemblyFile.out.contigs, wAssemblyFile.out.processed_reads)

    wUnmappedReadsList(wAssemblyFile.out.processed_reads, wBinning.out.bins)
    wFragmentRecruitmentList(wUnmappedReadsList.out.unmappedReads, Channel.fromPath(params?.steps?.fragmentRecruitment?.frhit?.genomes))

    wMagAttributesList(wBinning.out.bins)
    mapJoin(wMagAttributes.out.checkm, wBinning.out.bins_stats, "BIN_ID") | set { binsStats  }

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
