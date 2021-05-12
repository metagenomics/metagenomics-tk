nextflow.enable.dsl=2

params.out="out"
//params.input="/raid/simba/representatives_ids_dataset.tsv"
params.input="/raid/simba/bins_path_id_checkm_fixed_n50_fakedcoverage_final_checkm_filtered_id_dataset.tsv"
params.smetana_detailed=false

process carve {

    label 'tiny'
    tag "$sample $id"
    time '10h'

    errorStrategy 'ignore'
    publishDir "${params.out}/${sample}/carveme"
    input:
      tuple val(sample), val(id), path(mag_faa)
    output:
      tuple val("${sample}"), val("${id}"), path("${id}.xml"), emit: model
    shell:
    '''
    carve !{mag_faa} -o !{id}.xml
    '''
}

process memote {
    label 'tiny'

    errorStrategy 'ignore'

    tag "$sample $id"

    publishDir "${params.out}/${sample}/memote"

    input:
      tuple val(sample), val(id), path(model)

    output:
      tuple val("${sample}"), val("${id}"), path("${sample}_${id}_report.json.gz"), emit: report_json
      tuple val("${sample}"), val("${id}"), path("${sample}_${id}_report.html"), emit: report_html
      tuple val("${sample}"), val("${id}"), path("${sample}_${id}_metrics.tsv"), emit: report_tsv

    shell:
    '''
    memote run --solver cplex --filename !{sample}_!{id}_report.json.gz !{model}
    memote report snapshot --solver cplex --filename !{sample}_!{id}_report.html !{model}
    values=$(zcat !{sample}_!{id}_report.json.gz  | jq ' [ .tests.test_stoichiometric_consistency.duration, 
               .tests.test_reaction_mass_balance.metric, 
               .tests.test_reaction_charge_balance.metric, 
               .tests.test_find_disconnected.metric, 
               .tests.test_find_reactions_unbounded_flux_default_condition.metric ] | @tsv ')
    title=$(zcat !{sample}_!{id}_report.json.gz  | jq ' [ .tests.test_stoichiometric_consistency.title, 
               .tests.test_reaction_mass_balance.title, 
               .tests.test_reaction_charge_balance.title, 
               .tests.test_find_disconnected.title, 
               .tests.test_find_reactions_unbounded_flux_default_condition.title ] | @tsv ')
    echo -e "!{id}\t$title" > !{sample}_!{id}_metrics.tsv
    echo -e "!{id}\t$values" >> !{sample}_!{id}_metrics.tsv

    '''
}

process smetana_detailed {

    label 'large'

    tag "$sample"

    errorStrategy 'ignore'

    publishDir "${params.out}/${sample}/smetana/detailed/"

    when params.smetana_detailed

    input:
      tuple val(sample), path(xmls) 

    output:
      path("${sample}_detailed.tsv")

    shell:
    '''
    smetana !{xmls} -d -o !{sample}
    '''
}

process smetana_global {

    label 'large'

    tag "$sample"

    errorStrategy 'ignore'

    publishDir "${params.out}/${sample}/smetana/global/"

    input:
      tuple val(sample), path(xmls) 

    output:
      path("${sample}_global.tsv")

    shell:
    '''
    smetana !{xmls} -o !{sample}
    '''
}

process analyse {

    label 'tiny'
    tag "$sample $id"

    publishDir "${params.out}/${sample}/metabolites/"

    input:
      tuple val(sample), val(id), path(mag_json)  

    output:
      tuple val("${sample}"), val("${id}"), path("*.tsv"), emit: model

    shell:
    '''
    getProducts.sh !{mag_json} > !{id}_products.tsv
    getSubstrats.sh !{mag_json} > !{id}_substrats.tsv
    getReactions.sh !{mag_json} > !{id}_reactions.tsv
    '''
}


process build_json {
    tag "$sample $id"
    label 'tiny'
    publishDir "${params.out}/${sample}/json_output"
    errorStrategy 'retry'
    input:
      tuple val(sample), val(id), path(mag_xml)
    output:
      tuple val("${sample}"), val("${id}"), path("${id}.json"), emit: model
    shell:
    '''
    sbml_to_json.py !{mag_xml} !{id}.json cache
    '''
}


process prodigal {
    tag "$sample $id"
    label 'tiny'
    publishDir "${params.out}/${sample}/prodigal"
    input:
      tuple val(sample), val(id), path(mag)
    output:
      tuple val("${sample}"), val("${id}"), path("*.faa"), emit: protein
    shell:
    '''
    prodigal -i !{mag} -a  !{mag}.faa
    '''
}

workflow analyse_metabolites_file {
   take:
     bins
   main:
    Channel.from(file(bins)) | analyse_metabolites

   emit:
     carve_xml = carve.out.model
}

workflow analyse_metabolites {
   take:
     bins
   main:
//       | filter({ it.COMPLETENESS.toFloat() > 49 }) \
//       | filter({ it.CONTAMINATION.toFloat() < 5 })
       bins | map { it -> [it.SAMPLE, it.BIN_ID, it.PATH]}
       | set{binsChannel}

     binsChannel | prodigal | carve | build_json | analyse 
     carve.out | memote  
     carve.out | groupTuple(by:0) | map{ it -> it[0,2]} | set { model_group } 
     model_group | smetana_detailed 
     model_group | smetana_global

   emit:
     carve_xml = carve.out.model
}


workflow {
   analyse_metabolites(params.input) 
}
