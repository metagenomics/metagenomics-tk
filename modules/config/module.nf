nextflow.enable.dsl=2
import groovy.json.*
import org.yaml.snakeyaml.Yaml


def getOutput(SAMPLE, RUNID, TOOL, filename){
    return SAMPLE + '/' + RUNID + '/' + params.modules.config.name + '/' +
         params.modules.config.version.major + "." + 
         params.modules.config.version.minor + "." + 
         params.modules.config.version.patch +
         '/' + TOOL + '/' + filename
}


process pConfigUpload {

  label 'tiny'

  tag "$sample"

  publishDir params.output, mode: "${params.publishDirMode}", saveAs: { filename -> getOutput("${sample}", params.runid, "parameters", filename) }

  input:
  tuple val(sample), val(config), val(manifest)

  output:
  tuple val("${sample}"), file("*.yml")
  tuple file(".command.sh"), file(".command.out"), file(".command.err"), file(".command.log")

  shell:
  Yaml yaml = new Yaml()
  parameterName = config.steps.keySet().join(".")
  configStr = yaml.dump(config)
  workflowVersion = manifest.version
  workflowName = manifest.name
  '''
  echo "!{configStr}" > params_!{parameterName}.yml
  echo "version: !{workflowVersion}" >> manifest_!{parameterName}.yml
  echo "name: !{workflowName}" >> manifest_!{parameterName}.yml
  '''
}


/*
 * Takes a channel that must contain sample identifiers as values.
 *
 */
workflow wSaveSettingsList {
     take:
       samples
     main:
       samples | combine(Channel.from(params)) \
           | combine(Channel.from(workflow.manifest)) | pConfigUpload
}
