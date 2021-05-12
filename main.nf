nextflow.enable.dsl=2

include { analyse_metabolites } from './metabolomics/module'
include { bwa } from './read_mapping/bwa/module'
include { assembly_binning_input } from './assembly_binning/module'
include { dereplicate_file; dereplicate_list } from './dereplication/pasolli/module'

//readmapping
params.mapping_samples
params.list_of_representatives

workflow {
    analyse_metabolites(params.input)
}

workflow run_dereplication {
    dereplicate_file(Channel.from(file(params.input)))
}

workflow run_bwa {
    bwa(Channel.from('1'), Channel.from(params.input), Channel.fromPath(params.mapping_samples),Channel.fromPath(params.list_of_representatives))
}

workflow run_assembly_binning {
    assembly_binning_input(Channel.fromPath(params.input))
}

workflow run_pipeline {
    assembly_binning_input(Channel.fromPath(params.input))

    assembly_binning_input.out.processed_reads \
        | map { it -> "${it[0]}\t${it[1]}" } \
        | collectFile(seed: "SAMPLE\tREADS", name: 'test.txt', newLine: true) \
        | set { samples }

    dereplicate_list(assembly_binning_input.out.bins)
    representatives_list = dereplicate_list.out
    representatives_list | splitCsv(sep: '\t') \
       |  map { it -> file(it[0]) } | collectFile(){ item -> [ "representatives.fasta", item.text ] } | set { representatives_fasta }  

    bwa(Channel.from('1'), representatives_fasta, samples, representatives_list)

    assembly_binning_input.out.bins | map { it -> "${it.BIN_ID}\t${it.PATH}\t${it.SAMPLE}" } \
       | collectFile(newLine: true, seed: "BIN_ID\tPATH\tSAMPLE") | set { bin_attributes}

    assembly_binning_input.out.bins | analyse_metabolites
}
