nextflow.enable.dsl=2

def getOutput(SAMPLE, RUNID, TOOL, filename){
    return SAMPLE + '/' + 'RUNID' + '/' + params.modules.metabolomics.name + '/' +
          params.modules.metabolomics.version.major + "." + 
          params.modules.metabolomics.version.minor + "." +
          params.modules.metabolomics.version.patch + 
          '/' + TOOL + '/' + filename
}

process pCarveMe {

    label 'tiny'
    tag "$sample $id"
    time '10h'
    when params.steps.containsKey("metabolomics")

    errorStrategy 'ignore'
    publishDir params.output, saveAs: { filename -> getOutput("${sample}", params.runid, "carveme", filename) }

    input:
      tuple val(sample), val(id), path(mag_faa)
    output:
      tuple val("${sample}"), val("${id}"), path("${id}.xml.gz"), emit: model
    shell:
    '''
    carve !{mag_faa} -o !{id}.xml.gz
    '''
}

process pMemote {
    label 'tiny'

    errorStrategy 'ignore'

    tag "$sample $id"

    publishDir params.output, saveAs: { filename -> getOutput("${sample}", params.runid, "memote", filename) }

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

process pSmetanaDetailed {

    label 'large'

    tag "$sample"

    errorStrategy 'ignore'

    publishDir params.output, saveAs: { filename -> getOutput("${sample}", params.runid , "smetana/detailed/", filename) }

    when params?.steps?.metabolomics?.smetana?.contains("detailed")

    input:
      tuple val(sample), path(xmls) 

    output:
      path("${sample}_detailed.tsv")

    shell:
    '''
    smetana !{xmls} -d -o !{sample}
    '''
}

process pSmetanaGlobal {

    label 'large'
    tag "$sample"
    errorStrategy 'ignore'
    when params?.steps?.metabolomics?.smetana?.contains("global")

    publishDir params.output, saveAs: { filename -> getOutput("${sample}", params.runid, "smetana/global/", filename) }

    input:
      tuple val(sample), path(xmls) 

    output:
      path("${sample}_global.tsv")

    shell:
    '''
    smetana !{xmls} -o !{sample}
    '''
}

process pAnalyse {

    label 'tiny'
    tag "$sample $id"

    publishDir params.output, saveAs: { filename -> getOutput("${sample}", params.runid, "gsmmTsv", filename) }

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


process pBuildJson {
    tag "$sample $id"
    label 'tiny'
    publishDir params.output, saveAs: { filename -> getOutput("${sample}", params.runid, "gsmmJson", filename) }
    errorStrategy 'retry'
    input:
      tuple val(sample), val(id), path(mag_xml)
    output:
      tuple val("${sample}"), val("${id}"), path("${id}.json.gz"), emit: model
    shell:
    '''
    sbml_to_json.py !{mag_xml} !{id}.json
    gzip --best !{id}.json
    '''
}


process pProdigal {
    tag "$sample $id"
    label 'tiny'
    when params?.steps?.metabolomics !== null
    publishDir params.output, saveAs: { filename -> getOutput("${sample}", params.runid, "prodigal", filename) }
    input:
      tuple val(sample), val(id), path(mag)
    output:
      tuple val("${sample}"), val("${id}"), path("*.faa.gz"), emit: protein
    shell:
    '''
    cat !{mag} | prodigal -a !{mag}.faa
    gzip --best !{mag}.faa
    '''
}

workflow wAnalyseMetabolitesFile {
   take:
     bins
   main:
    Channel.from(file(bins)) | analyse_metabolites
   emit:
     carve_xml = carve.out.model
}

workflow wAnalyseMetabolites {
   take:
     bins
   main:
//       | filter({ it.COMPLETENESS.toFloat() > 49 }) \
//       | filter({ it.CONTAMINATION.toFloat() < 5 })
       bins | map { it -> [it.SAMPLE, it.BIN_ID, it.PATH]}
       | set{binsChannel}

     binsChannel | pProdigal | pCarveMe | pBuildJson | pAnalyse 
     pCarveMe.out | pMemote  
     pCarveMe.out | groupTuple(by:0) | map{ it -> it[0,2]} | set { model_group } 
     model_group | pSmetanaDetailed 
     model_group | pSmetanaGlobal

   emit:
     carve_xml = pCarveMe.out.model
}
