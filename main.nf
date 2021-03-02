nextflow.enable.dsl=2

include { analyse_metabolites } from './metabolomics/module'

workflow {
    analyse_metabolites(params.input)
}
