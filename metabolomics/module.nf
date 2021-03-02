nextflow.enable.dsl=2

params.out="out"
//params.input="/raid/simba/representatives_ids_dataset.tsv"
params.input="/raid/simba/bins_path_id_checkm_fixed_n50_fakedcoverage_final_checkm_filtered_id_dataset.tsv"

process carve {
    cpus 2
    publishDir params.out + "/carveme"
    conda '/vol/spool/CAMI/CAMI_MOUSEGUT/test/almeida/conda/envs/carveme_env'
    input:
      tuple path(mag_faa), val(id)
    output:
      tuple path("${id}.xml"), val("${id}"), emit: model
    shell:
    '''
    carve !{mag_faa} -o !{id}.xml
    '''
}


process smetana_detailed {
    cpus 14
    errorStrategy 'ignore'
    publishDir params.out + "/smetana/detailed"
    conda '/vol/spool/CAMI/CAMI_MOUSEGUT/test/almeida/conda/envs/carveme_env'
    input:
      tuple path(xmls), val(dataset_id)
    output:
      path("${dataset_id}_detailed.tsv")
    shell:
    '''
    smetana !{xmls} -d -o !{dataset_id}
    '''
}

process smetana_global {
    cpus 14
    errorStrategy 'ignore'
    publishDir params.out + "/smetana/global"
    conda '/vol/spool/CAMI/CAMI_MOUSEGUT/test/almeida/conda/envs/carveme_env'
    input:
      tuple path(xmls), val(dataset_id)
    output:
      path("${dataset_id}_global.tsv")
    shell:
    '''
    smetana !{xmls} -o !{dataset_id}
    '''
}

process analyse {
    conda 'jq'
    cpus 2
    publishDir params.out + "/metabolites"
    conda '/vol/spool/CAMI/CAMI_MOUSEGUT/test/almeida/conda/envs/carveme_env'
    input:
      tuple path(mag_json), val(id)
    output:
      tuple path("*.tsv"), val("${id}"), emit: model
    shell:
    '''
    getProducts.sh !{mag_json} > !{id}_products.tsv
    getSubstrats.sh !{mag_json} > !{id}_substrats.tsv
    getReactions.sh !{mag_json} > !{id}_reactions.tsv
    '''
}


process build_json {
    conda 'bioconda::cobra python=3.8'
    cpus 2
    publishDir params.out + "/json_output"
    input:
      tuple path(mag_xml), val(id)
    output:
      tuple path("${id}.json"), val("${id}"), emit: model
    shell:
    '''
    sbml_to_json.py !{mag_xml} !{id}.json
    '''
}


process prodigal {
    cpus 1
    publishDir params.out + "/prodigal"
    conda 'prodigal'
    input:
      tuple val(id), path(mag)
    output:
      tuple path("*.faa"), val("${id}"), emit: protein
    shell:
    '''
    prodigal -i !{mag} -a  !{mag}.faa
    '''
}


workflow analyse_metabolites {
   take:
     bins
   main:
    Channel.from(file(bins)) | splitCsv(sep: '\t', header: true) \
       | filter({ it.COMPLETENESS.toFloat() > 49 }) \
       | filter({ it.CONTAMINATION.toFloat() < 5 })
       | map { it -> [it.DATASET, it.BIN, it.PATH]}
       | set{binsChannel}

//     Channel.from(file(bins)) | splitCsv(sep: "\t") | set{binsChannel}  
     binsChannel | map { it -> [it[1], it[2]]} | prodigal | carve | build_json | analyse 
     carve.out | join(binsChannel, by:1) | groupTuple(by:2) | map{ it -> it[1,2]} | (smetana_detailed & smetana_global) 
   emit:
     carve_xml = carve.out.model
}


workflow {
   analyse_metabolites(params.input) 
}
