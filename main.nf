nextflow.enable.dsl=2


include { wAssemblyFile } from './assembly/module'
include { wBinning } from './binning/module.nf'
include { wMagAttributes } from './magAttributes/module.nf'
include { wDereplicateFile; wDereplicateList } from './dereplication/pasolli/module'
include { wReadMappingBwa } from './readMapping/bwa/module'
include { wAnalyseMetabolites } from './metabolomics/module'


def mapJoin(channel_a, channel_b, key){
    channel_a \
        | map{ it -> [it[key], it] } \
        | cross(channel_b | map{it -> [it[key], it]}) \
        | map { it[0][1] + it[1][1] }
}


workflow {
//   analyse_metabolites(params.input)
}


workflow run_dereplication {
//   dereplicate_file(Channel.from(file(params.input)))
}


workflow run_bwa {
//    bwa(Channel.from('1'), Channel.from(params.input), Channel.fromPath(params.mapping_samples),Channel.fromPath(params.list_of_representatives))
}


workflow mags_generation {
//    mags_generation_file(Channel.fromPath(params.input))
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
    file(params.tempdir).mkdirs()

    wAssemblyFile(Channel.fromPath(params.input))

    wAssemblyFile.out.processed_reads \
        | map { it -> "${it[0]}\t${it[1]}" } \
        | collectFile(seed: "SAMPLE\tREADS", name: 'test.txt', newLine: true) \
        | set { samples }

    wBinning(wAssemblyFile.out.contigs, wAssemblyFile.out.processed_reads)

    wMagAttributes(wBinning.out.bins, wBinning.out.mapping)
    mapJoin(wMagAttributes.out.bins_info, wBinning.out.bins_stats, "BIN_ID") | set { binsStats  }

    wDereplicateList(binsStats)
    representativesList = wDereplicateList.out

    REPRESENTATIVES_PATH_IDX = 0
    representativesList | splitCsv(sep: '\t') \
       |  map { it -> file(it[REPRESENTATIVES_PATH_IDX]) } \
       | collectFile(tempDir: file(params.tempdir)){ item -> [ "representatives.fasta", item.text ] } \
       | set { representativesFasta }  

    wReadMappingBwa(Channel.from('1'), representativesFasta, samples, representativesList)

    binsStats | wAnalyseMetabolites
}
