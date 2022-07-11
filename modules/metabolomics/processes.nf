def getOutput(SAMPLE, RUNID, TOOL, filename){
    return SAMPLE + '/' + RUNID + '/' + params.modules.metabolomics.name + '/' +
          params.modules.metabolomics.version.major + "." + 
          params.modules.metabolomics.version.minor + "." +
          params.modules.metabolomics.version.patch + 
          '/' + TOOL + '/' + filename
}

process pCarveMe {

    label 'tiny'

    tag "Sample: $sample, Bin: $id"

    container "${params.carveme_image}"

    containerOptions " --user 0:0 "

    when params.steps.containsKey("metabolomics") && params.steps.metabolomics.containsKey("carveme")

    beforeScript  params?.steps.containsKey("metabolomics") \
	? Utils.getBeforeScript(params?.steps?.metabolomics?.beforeProcessScript.trim(), params.carveme_image) \
	: ""
  
    publishDir params.output, saveAs: { filename -> getOutput("${sample}", params.runid, "carveme", filename) }, \
        pattern: "{**.xml}"

    input:
      tuple val(sample), val(id), path(mag)
      val(mode)

    output:
      tuple val("${sample}"), val("${id}"), path("${id}.xml"), optional: true, emit: model
      tuple val("${id}"), val("${output}"), val(params.LOG_LEVELS.INFO), file(".command.sh"), \
	file(".command.out"), file(".command.err"), file(".command.log"), emit: logs

    shell:
    output = getOutput("${sample}", params.runid, "carveme", "")

    if(mode == "proteins")
       '''
       carve !{mag} -o !{id}.model.xml !{params.steps.metabolomics.carveme.additionalParams}
       '''
    else if(mode == "genome")
       '''
       carve --dna !{mag} -o !{id}.model.xml !{params.steps.metabolomics.carveme.additionalParams}
       '''
     else
        error "Invalid mode: ${mode}"
}
