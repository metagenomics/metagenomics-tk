nextflow.enable.dsl=2

include { analyse_metabolites } from './metabolomics/module'
include { bwa } from './read_mapping/bwa/module'

params.mapping_samples

workflow {
    analyse_metabolites(params.input)
}

workflow bwa {
    bwa(Channel.from('1'), Channel.from(params.input), Channel.fromPath(params.mapping_samples))
}
