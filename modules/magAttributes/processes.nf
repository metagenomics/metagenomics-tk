def getOutput(SAMPLE, RUNID, TOOL, filename){
    return SAMPLE + '/' + RUNID + '/' + params.modules.magAttributes.name + '/' + 
         params.modules.magAttributes.version.major + "." +  
         params.modules.magAttributes.version.minor + "." +
         params.modules.magAttributes.version.patch +
         '/' + TOOL + '/' + filename
}

/**
*
* The pCollectFile process is designed to collect and concatenate files of checkm and gtdb.
* It uses the csvtk tool to concatenate all the .tsv files.
*
**/
process pCollectFile {

    memory { Utils.getMemoryResources(params.resources.small, "${sample}", task.attempt, params.resources) }

    cpus { Utils.getCPUsResources(params.resources.small, "${sample}", task.attempt, params.resources) }

    tag "Sample: $sample"

    container "${params.ubuntu_image}"

    publishDir params.output, mode: "${params.publishDirMode}", saveAs: { filename -> getOutput("${sample}",params.runid ,"${module}", filename) }, \
      pattern: "{**.tsv}"

    input:
    val(module)
    val(tool)
    val(type)
    val(name)
    tuple val(sample), path(gtdbSummaries)

    output:
    tuple val("${sample}"), path("${sample}_${tool}_${type}${name}")

    shell:
    '''
    csvtk -t -T concat !{gtdbSummaries} > !{sample}_!{tool}_!{type}!{name}
    '''
}


workflow _wCollectFile {
  take:
    type
    inputFile
    module
    tool
    nameSuffix
  main:
    inputFile | map { sample, path, size -> tuple( groupKey(sample, size), path ) } \
	| groupTuple(remainder: true) | set {collectedCombinedFiles}
    pCollectFile(module, tool, type, nameSuffix, collectedCombinedFiles)
    pCollectFile.out | set { combinedListFiles }
  emit:
    combinedListFiles
}
