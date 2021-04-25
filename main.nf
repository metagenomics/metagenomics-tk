nextflow.enable.dsl=2

include { analyse_metabolites } from './metabolomics/module'
include { bwa } from './read_mapping/bwa/module'
include { assembly_binning_input } from './assembly_binning/module'

//readmapping
params.mapping_samples
params.list_of_representatives

workflow {
    analyse_metabolites(params.input)
}

workflow run_bwa {
    bwa(Channel.from('1'), Channel.from(params.input), Channel.fromPath(params.mapping_samples),Channel.fromPath(params.list_of_representatives))
}

workflow run_assembly_binning {
    assembly_binning_input(Channel.fromPath(params.input))
}
