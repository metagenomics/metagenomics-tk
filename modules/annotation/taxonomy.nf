include { pMMseqs2_taxonomy } from './mmseqs2'
include { pDumpLogs } from '../utils/processes'
include { runMMSeqsTaxonomy } from '../utils/methods'

workflow _wCreateTaxonomy {
  take:
    sourceChannel
    contigs
    proteins
    fastaCounter
    contig2GeneMapping
  main:

    SAMPLE_IDX=0
    PATH_IDX=2

    // MMseqs2 taxonomy and Metabuli both write to taxonomy/emgb and feed the EMGB export, so only one may be configured.
    if (params?.steps?.annotation?.containsKey("mmseqs2_taxonomy") && params?.steps?.annotation?.containsKey("metabuli")) {
        error "Configure either annotation.mmseqs2_taxonomy or annotation.metabuli, not both."
    }

    selectedMMSeqsTaxDBs = (params?.steps?.annotation?.mmseqs2_taxonomy ?: [:])
        .findAll { it.key != "runOnMAGs" }
        .collect {
            [
                it.key,
                it.value?.params ?: "",
                it.value?.ramMode ? "true" : "false",
                it.value?.initialMaxSensitivity ?: "",
                it.value?.database?.extractedDBPath ?: "",
                it.value.database?.download?.source ?: "",
                it.value.database?.download?.md5sum ?: "",
                it.value.database?.download?.s5cmd?.params ?: "",
            ]
        }

    proteins \
        | map { [it[SAMPLE_IDX], it[PATH_IDX]] }
        | combine(fastaCounter, by: SAMPLE_IDX)
        | map { sample, path, size -> tuple(groupKey(sample, size), path) }
        | groupTuple()
        | combine(Channel.from(selectedMMSeqsTaxDBs))
        | set { combinedMMseqsTax }


    selectedMetabuliTaxDBs = (params?.steps?.annotation?.metabuli ?: [:])
        .findAll { it.key != "runOnMAGs" }
        .collect {
            [
                it.key,
                it.value?.params ?: "",
                it.value?.database?.extractedDBPath ?: "",
                it.value.database?.download?.source ?: "",
                it.value.database?.download?.md5sum ?: "",
                it.value.database?.download?.s5cmd?.params ?: "",
            ]
        }

    proteins \
        | map { [it[SAMPLE_IDX], it[PATH_IDX]] }
        | combine(fastaCounter, by: SAMPLE_IDX)
        | map { sample, path, size -> tuple(groupKey(sample, size), path) }
        | groupTuple()
        | combine(Channel.from(selectedMetabuliTaxDBs))
        | set { combinedMetabuliTax }

    contigs | map { sample, binId, paths, kingdom -> [sample, binId, paths]} 
        | combine(Channel.from(selectedMetabuliTaxDBs))
        | set { combinedMetabuliTaxInput }

    pMMseqs2_taxonomy(sourceChannel, combinedMMseqsTax)

    pMetabuliTaxonomy(sourceChannel, combinedMetabuliTaxInput)

    pMetabuliTaxonomy.out.classifications
        | map { sample, binType, binId, dbType, path -> [sample, binType, dbType, path] }
        | combine(fastaCounter, by: SAMPLE_IDX)
        | map { sample, binType, dbType, path, size -> tuple(groupKey([sample, binType, dbType], size), path) }
        | groupTuple()
        | map { key, paths -> [*key.getGroupTarget(), paths] }
        | combine(contig2GeneMapping, by: SAMPLE_IDX)
        | pMetabuliToEMGB

    pMMseqs2_taxonomy.out.logs 
        | mix(pMetabuliTaxonomy.out.logs)
        | mix(pMetabuliToEMGB.out.logs)
        | pDumpLogs

    pMMseqs2_taxonomy.out.taxonomy 
        | mix(pMetabuliToEMGB.out.taxonomy)
        | set {taxonomy}

    pMMseqs2_taxonomy.out.krakenStyleTaxonomy
        | mix(pMetabuliTaxonomy.out.krakenReport)
        | set { krakenReport }

 emit:
    kronaHtml = kronaHtml
    mmseqs2_krakenTaxonomy = pMMseqs2_taxonomy.out.krakenReport
    taxonomy = taxonomy
}

/**
*
* Metabuli classifies contigs or MAGs taxonomically by comparing their DNA sequences using metamers.
*
**/
process pMetabuliTaxonomy {
   
      container "${params.metabuli_image}"

      containerOptions Utils.constructParametersObject("metabuli", params, params.apptainer) + (params.apptainer ? "" : " --entrypoint='' " + Utils.getDockerNetwork())

      tag "Sample: $sample, BinID: $binId, Database: $dbType"

      memory { Utils.getMemoryResources(params.resources.large, "${sample}", task.attempt, params.resources) }

      cpus { Utils.getCPUsResources(params.resources.large, "${sample}", task.attempt, params.resources) }

      secret { "${S3_TAX_DB_ACCESS}"!="" ? ["S3_TAX_${dbType}_ACCESS", "S3_TAX_${dbType}_SECRET"] : [] } 

      publishDir params.output, mode: "${params.publishDirMode}", saveAs: { filename -> Output.getOutput("${sample}", params.runid, "taxonomy/metabuli/${dbType}", params.modules.annotation, new File(filename).getName()) }, \
         pattern: "out/*{_classifications.tsv,_report.tsv,_krona.html}"
 
      when:
      runMMSeqsTaxonomy(params?.steps.containsKey("annotation") && params?.steps.annotation.containsKey("metabuli"), \
	   binType.equals("binned"), \
	   binType.equals("plasmid"), \
           params?.steps.containsKey("annotation") && params?.steps.annotation.containsKey("metabuli") && params?.steps.annotation.metabuli.runOnMAGs)

   input:
      val(binType)
      tuple val(sample), val(binId), file(contigs), val(dbType), val(parameters), val(EXTRACTED_DB), val(DOWNLOAD_LINK), val(MD5SUM), val(S5CMD_PARAMS)

   output:
      tuple val("${sample}"), val("${binType}"), val("${binId}"), val("${dbType}"), path("out/*_classifications.tsv"), optional:true, emit: classifications
      tuple val("${sample}"), val("${dbType}"), path("out/*_report.tsv"), optional:true, emit: krakenReport 
      tuple val("${sample}"), val("${dbType}"), path("out/*_krona.html"), optional:true, emit: kronaHtml
      tuple val("${sample}_${binType}"), val("${output}"), val(params.LOG_LEVELS.INFO), file(".command.sh"), \
        file(".command.out"), file(".command.err"), file(".command.log"), emit: logs

   shell:
   output = Output.getOutput("${sample}", params.runid, "taxonomy/metabuli/${dbType}", params.modules.annotation, "")
   S3_TAX_DB_ACCESS=params.steps?.annotation?.metabuli?."${dbType}"?.database?.download?.s5cmd && S5CMD_PARAMS.indexOf("--no-sign-request") == -1 ? "\$S3_TAX_${dbType}_ACCESS" : ""
   S3_TAX_DB_SECRET=params.steps?.annotation?.metabuli?."${dbType}"?.database?.download?.s5cmd && S5CMD_PARAMS.indexOf("--no-sign-request") == -1 ? "\$S3_TAX_${dbType}_SECRET" : ""
   template("metabuli.sh")
}


/**
*
* Merges the Metabuli classifications of all bins of a sample and converts them into the MMseqs2 taxonomy TSV format,
* so that downstream steps (e.g. the EMGB export) can use both tools interchangeably.
*
**/
process pMetabuliToEMGB {

      container "${params.metabuli_image}"

      containerOptions Utils.constructParametersObject("metabuli", params, params.apptainer) + (params.apptainer ? "" : " --entrypoint='' " + Utils.getDockerNetwork())

      tag "Sample: $sample, Database: $dbType"

      memory { Utils.getMemoryResources(params.resources.tiny, "${sample}", task.attempt, params.resources) }

      cpus { Utils.getCPUsResources(params.resources.tiny, "${sample}", task.attempt, params.resources) }

      publishDir params.output, mode: "${params.publishDirMode}", \
         saveAs: { filename -> Output.getOutput("${sample}", params.runid, "taxonomy/emgb/${dbType}", params.modules.annotation, new File(filename).getName()) }, \
         pattern: "*.taxonomy.tsv"

   input:
      tuple val(sample), val(binType), val(dbType), path(classifications), path(contig2GeneMapping)

   output:
      tuple val("${sample}"), val("${dbType}"), path("*.taxonomy.tsv"), emit: taxonomy
      tuple val("${sample}_${binType}"), val("${output}"), val(params.LOG_LEVELS.INFO), file(".command.sh"), \
        file(".command.out"), file(".command.err"), file(".command.log"), emit: logs

   shell:
   output = Output.getOutput("${sample}", params.runid, "taxonomy/emgb/${dbType}", params.modules.annotation, "")
   template("metabuliToEMGB.sh")
}
