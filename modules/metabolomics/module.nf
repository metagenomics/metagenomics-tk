include { pDumpLogs } from '../utils/processes'
include { wSaveSettingsList } from '../config/module'

include { pCarveMe } from './processes'

def getOutput(SAMPLE, RUNID, TOOL, filename){
    return SAMPLE + '/' + RUNID + '/' + params.modules.metabolomics.name + '/' +
          params.modules.metabolomics.version.major + "." + 
          params.modules.metabolomics.version.minor + "." +
          params.modules.metabolomics.version.patch + 
          '/' + TOOL + '/' + filename
}

process pGapSeq {

    label 'small'

    tag "Sample: $sample, Bin: $id"

    beforeScript params?.steps.containsKey("metabolomics") \
	? Utils.getBeforeScript(params?.steps?.metabolomics?.beforeProcessScript.trim(), params.gapseq_image) \
	: ""

    when params.steps.containsKey("metabolomics") \
	&& params.steps.metabolomics.containsKey("gapseq") \
	&& !params.steps.metabolomics.containsKey("carveme")

    container "${params.gapseq_image}"

    publishDir params.output, mode: "${params.publishDirMode}", \
      saveAs: { filename -> getOutput("${sample}", params.runid, "gapseq", filename) }, \
      pattern: "{**.xml,**.tbl,**.RDS,**.csv}"

    input:
      tuple val(sample), val(id), path(mag)

    output:
      tuple val("${sample}"), val("${id}"), path("*-draft.xml"), optional: true, emit: draft
      tuple val("${sample}"), val("${id}"), path("${id}.model.xml"), optional: true, emit: model
      tuple val("${sample}"), val("${id}"), path("*.tbl"), optional: true, emit: tables
      tuple val("${sample}"), val("${id}"), path("*.RDS"), optional: true, emit: rds
      tuple val("${sample}"), val("${id}"), path("*.csv"), optional: true, emit: csv
      tuple val("${id}"), val("${output}"), val(params.LOG_LEVELS.INFO), file(".command.sh"), \
	file(".command.out"), file(".command.err"), file(".command.log"), emit: logs

    shell:
    output = getOutput("${sample}", params.runid, "gapseq", "")
    '''
    gapseq doall !{mag} !{params.steps.metabolomics.gapseq.additionalParams}
    mv $(ls -1 *.xml | grep -v -- "-draft.xml") !{id}.model.xml
    '''
}


process pMemote {

    label 'highmemMedium'

    tag "Sample: $sample, Bin: $id"

    container "${params.memote_image}"

    beforeScript params?.steps.containsKey("metabolomics") \
	? Utils.getBeforeScript(params?.steps?.metabolomics?.beforeProcessScript.trim(), params.memote_image) \
	: ""

    when params.steps.containsKey("metabolomics") && params.steps.metabolomics.containsKey("memote")

    publishDir params.output, mode: "${params.publishDirMode}", saveAs: { filename -> getOutput("${sample}", params.runid, "memote", filename) }, \
      pattern: "{**.json.gz,**.html,**.tsv}"

    input:
      tuple val(sample), val(id), path(model)

    output:
      tuple val("${sample}"), val("${id}"), path("${sample}_${id}_report.json.gz"), optional: true, emit: reportJson
      tuple val("${sample}"), val("${id}"), path("${sample}_${id}_report.html"), optional: true, emit: reportHtml
      tuple val("${sample}"), val("${id}"), path("${sample}_${id}_metrics.tsv"), optional: true, emit: reportTsv
      tuple val("${id}"), val("${output}"), val(params.LOG_LEVELS.INFO), file(".command.sh"), \
	file(".command.out"), file(".command.err"), file(".command.log"), emit: logs

    shell:
    output = getOutput("${sample}", params.runid, "memote", "")
    template "memote.sh"
}

process pSmetanaDetailed {

    label 'highmemLarge'

    tag "Sample: $sample"

    container "${params.smetana_image}"

    beforeScript params?.steps.containsKey("metabolomics") \
	? Utils.getBeforeScript(params?.steps?.metabolomics?.beforeProcessScript.trim(), params.smetana_image) \
	: ""

    publishDir params.output, mode: "${params.publishDirMode}", saveAs: { filename -> getOutput("${sample}", params.runid , "smetana/detailed/", filename) }, \
      pattern: "{**.tsv}"

    when params?.steps.containsKey("metabolomics") && params?.steps?.metabolomics?.containsKey("smetana") \
       && params?.steps?.metabolomics?.smetana.containsKey("detailed") \
       && params?.steps?.metabolomics?.smetana?.detailed

    input:
      tuple val(sample), path(xmls) 

    output:
      path("${sample}_detailed.tsv")
      tuple file(".command.sh"), file(".command.out"), file(".command.err"), file(".command.log")

    shell:
    template "smetanaDetailed.sh"
}


process pSmetanaGlobal {

    label 'highmemLarge'

    tag "Sample: $sample"

    beforeScript params?.steps.containsKey("metabolomics") \
	? Utils.getBeforeScript(params?.steps?.metabolomics?.beforeProcessScript.trim(), params.smetana_image) \
	: ""

    container "${params.smetana_image}"

    when params?.steps.containsKey("metabolomics") && params?.steps?.metabolomics?.containsKey("smetana") \
       && params?.steps?.metabolomics?.smetana.containsKey("detailed") \
       && params?.steps?.metabolomics?.smetana?.global

    publishDir params.output, mode: "${params.publishDirMode}", saveAs: { filename -> getOutput("${sample}", params.runid, "smetana/global/", filename) }, \
      pattern: "{**.tsv}"

    input:
      tuple val(sample), path(xmls) 

    output:
      path("${sample}_global.tsv")
      tuple file(".command.sh"), file(".command.out"), file(".command.err"), file(".command.log")

    shell:
    template "smetanaGlobal.sh"
}


process pAnalyse {

    label 'tiny'

    tag "Sample: $sample, Bin: $id"

    container "${params.ubuntu_image}"

    publishDir params.output, mode: "${params.publishDirMode}", \
       saveAs: { filename -> getOutput("${sample}", params.runid, "gsmmTsv", filename) }, \
       pattern: "{**.tsv}"

    input:
      tuple val(sample), val(id), path(magJson)  

    output:
      tuple val("${sample}"), val("${id}"), path("*.tsv"), emit: model
      tuple val("${id}"), val("${output}"), val(params.LOG_LEVELS.INFO), file(".command.sh"), \
	file(".command.out"), file(".command.err"), file(".command.log"), emit: logs

    shell:
    output = getOutput("${sample}", params.runid, "gsmmTsv", "")
    '''
    PRODUCTS=!{id}_products.tsv
    REACTIONS=!{id}_reactions.tsv
    SUBSTRATS=!{id}_substrats.tsv

    zcat -f !{magJson}  \
	| jq -r '.|.reactions[] | .metabolites | to_entries | .[] | select(.value>=0) | .key ' |  tr -d '"' | sed -r '/^\s*$/d'  >  ${PRODUCTS} 
    sed -i -e '1i BIN_ID\tPRODUCTS' -e "s/^/!{id}\t/" ${PRODUCTS}

    zcat -f !{magJson} \
	| jq -r '.|.reactions[] |.name ' | sed -r '/^\s*$/d'  | sort > ${REACTIONS}
    sed -i -e '1i BIN_ID\tREACTIONS' -e "s/^/!{id}\t/" ${REACTIONS}

    zcat -f !{magJson} \
	| jq -r '.|.reactions[] | .metabolites | to_entries | .[] | select(.value<0) | .key ' |  tr -d '"' | sed -r '/^\s*$/d' > ${SUBSTRATS}
    sed -i -e '1i BIN_ID\tSUBSTRATS' -e "s/^/!{id}\t/" ${SUBSTRATS}
    '''
}


process pBuildJson {

    tag "Sample: $sample, Bin: $id"

    label 'tiny'

    container "${params.cobra_image}"

    containerOptions " --user 0:0"

    publishDir params.output, mode: "${params.publishDirMode}", \
      saveAs: { filename -> getOutput("${sample}", params.runid, "gsmmJson", filename) }, \
      pattern: "{**.json}"

    input:
      tuple val(sample), val(id), path(magXml)

    output:
      tuple val("${sample}"), val("${id}"), path("${id}.json"), emit: model
      tuple val("${id}"), val("${output}"), val(params.LOG_LEVELS.INFO), file(".command.sh"), \
	file(".command.out"), file(".command.err"), file(".command.log"), emit: logs

    shell:
    output = getOutput("${sample}", params.runid, "gsmmJson", "")
    '''
    sbml_to_json.py !{magXml} !{id}.json
    '''
}


workflow wAnalyseMetabolitesFile {
   main:
     bins = Channel.empty()
     proteins = Channel.empty()
     type = Channel.empty()
     if(params.steps.metabolomics.input.containsKey("proteins")){
         Channel.from(file(params.steps.metabolomics.input.proteins)) \
		| splitCsv(sep: '\t', header: true) \
		| map {it -> [ it.SAMPLE, it.BIN_ID, it.PATH ]} \
                | set { proteins }
         Channel.value("proteins") | set { type }
     }
     if(params.steps.metabolomics.input.containsKey("bins")){
         Channel.from(file(params.steps.metabolomics.input.bins))
		| splitCsv(sep: '\t', header: true) \
		| map {it -> [ it.SAMPLE, it.BIN_ID, it.PATH ]} \
                | set { bins }
         Channel.value("genome") | set { type }
     }

     SAMPLE_IDX = 0
     wSaveSettingsList(bins | mix(proteins) |  map { it[SAMPLE_IDX] } \
	| unique | map { it -> it[SAMPLE_IDX] })

     _wAnalyseMetabolites(bins, proteins, type)
  emit:
     models = _wAnalyseMetabolites.out.models
}


/*
*  This entrypoint accepts a channel containing dictionaries as values.
*  The bins channel must contain the following values: 
*  [COMPLETENESS:, CONTAMINATION:, SAMPLE:, BIN_ID:, PATH: ]
*
*  The proteins channel must contain the following values: 
*  [COMPLETENESS:, CONTAMINATION:, SAMPLE:, BIN_ID:, PROTEINS: ]
*
*  where:
*     - COMPLETENESS and CONTAMINATION are values extracted from the checkm output.
*     - BIN_ID is a unique id for the input genome or protein file over all samples.
*     - SAMPLE is a unique identifier for the sample
*     - PROTEINS and PATH contain the path to the corresponding genome or protein file.
*/
workflow wAnalyseMetabolitesList {
   take:
     bins
     proteins
   main:
     // Filter by completeness and contamination

     filteredProteins = Channel.empty()
     filteredBins = Channel.empty()
     type = Channel.empty()

     if(params.steps.containsKey("metabolomics") && params.steps.metabolomics.containsKey("gapseq")){
       bins | filter({ it.COMPLETENESS.toFloat() > params.steps?.metabolomics?.filter?.minCompleteness }) \
          | filter({ it.CONTAMINATION.toFloat() < params.steps?.metabolomics?.filter?.maxContamination }) \
          | map { it -> [it.SAMPLE, it.BIN_ID, it.PATH]} \
          | set { filteredBins }

       Channel.value("genome") | set { type }
     }

     if(params.steps.containsKey("metabolomics") && params.steps.metabolomics.containsKey("carveme")){
       proteins | filter({ it.COMPLETENESS.toFloat() > params.steps?.metabolomics?.filter?.minCompleteness }) \
          | filter({ it.CONTAMINATION.toFloat() < params.steps?.metabolomics?.filter?.maxContamination }) \
          | map { it -> [it.SAMPLE, it.BIN_ID, it.PROTEINS]} \
          | set { filteredProteins }

       Channel.value("proteins") | set { type }
     }

     _wAnalyseMetabolites(filteredBins, filteredProteins, type)
   emit:
     models = _wAnalyseMetabolites.out.models
}


workflow _wAnalyseMetabolites {
     take:
        bins
        proteins
        type
     main:
        bins | pGapSeq

        // While GapSeq is only able to process genomes, carveme is also able to use
        // predicted proteins.
        pCarveMe(proteins | mix(bins), type)

        // build, validate and analyse all models
        pGapSeq.out.model | mix(pCarveMe.out.model) | set { model } 

        model | pBuildJson & pMemote
        pBuildJson.out.model | pAnalyse 

        // group models per sample
        SAMPLE_IDX = 0
        MODEL_IDX = 2
        model | groupTuple(by: SAMPLE_IDX) \
	 | map{ it -> it[SAMPLE_IDX, MODEL_IDX]} | set { modelGroup } 

        // compute possible interactions 
        modelGroup | pSmetanaDetailed & pSmetanaGlobal

        pGapSeq.out.logs | mix(pCarveMe.out.logs) \
            | mix(pMemote.out.logs, pAnalyse.out.logs, pBuildJson.out.logs) | pDumpLogs
     emit:
        models = model
}
