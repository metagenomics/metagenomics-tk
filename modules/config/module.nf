def getYaml() {
  return new org.yaml.snakeyaml.Yaml()
}


process pConfigUpload {

  label 'tiny'

  tag "${sample}"

  publishDir params.output, mode: "${params.publishDirMode}", saveAs: { filename -> Output.getOutput("${sample}", params.runid, "parameters", params.modules.config, filename) }

  input:
  tuple val(sample), val(config), val(manifest)

  output:
  tuple val("${sample}"), file("*.yml")
  tuple file(".command.sh"), file(".command.out"), file(".command.err"), file(".command.log")

  script:
  configStr = getYaml().dump(config)
  workflowVersion = manifest.version
  workflowName = manifest.name
  timestamp = new java.util.Date().format('YYYYMMdd-HHmmss-SSS')
  """
  echo "${configStr}" > params_${timestamp}.yml
  echo "version: ${workflowVersion}" >> manifest_${timestamp}.yml
  echo "name: ${workflowName}" >> manifest_${timestamp}.yml
  """
}


/*
 * Takes a channel that must contain sample identifiers as values.
 *
 */
workflow wSaveSettingsList {
  take:
  samples

  main:
  samples
    | combine(channel.from(params))
    | combine(channel.from(workflow.manifest))
    | pConfigUpload
}
