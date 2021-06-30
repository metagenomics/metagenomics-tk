nextflow.enable.dsl=2
import groovy.json.*
import org.yaml.snakeyaml.Yaml


MODULE="config"
VERSION="0.1.0"
def getOutput(SAMPLE, RUNID, TOOL, filename){
    return SAMPLE + '/' + RUNID + '/' + MODULE + '/' + VERSION + '/' + TOOL + '/' + filename
}


process pConfigUpload {

  label 'tiny'

  tag "$sample"

  publishDir params.output, saveAs: { filename -> getOutput("${sample}", params.runid, "parameters", filename) }

  input:
  tuple val(sample), val(config), val(manifest)

  output:
  tuple val("${sample}"), file("*.yml")

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
 * Takes a tab separated file of files containing sample identifiers as input and places meta-omics-toolkit specific files in sample directories.
 * Input tsv file must contain a column with SAMPLE as header.
 *
 */
workflow wSaveSettingsFile {
     take:
       samplesTable
     main:
          
        // Parse TSV file to get sample ids 
        samplesTable | splitCsv(sep: '\t', header: true) \
           | map { it -> it.SAMPLE} | set{ samples }  \

       samples | combine(Channel.from(params)) \
           | combine(Channel.from(workflow.manifest)) | pConfigUpload
}
