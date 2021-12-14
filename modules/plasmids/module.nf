nextflow.enable.dsl=2

include { pDumpLogs } from '../utils/processes'

def getOutput(SAMPLE, RUNID, TOOL, filename){
    return SAMPLE + '/' + RUNID + '/' + params.modules.plasmids.name + '/' +
          params.modules.plasmids.version.major + "." +
          params.modules.plasmids.version.minor + "." +
          params.modules.plasmids.version.patch +
          '/' + TOOL + '/' + filename
}


process pSCAPP {

    label 'medium'

    tag "$sample"

    publishDir params.output, saveAs: { filename -> getOutput("${sample}", params.runid, "SCAPP", filename) }

    when params.steps.containsKey("plasmid") && params.steps.plasmid?.containsKey("SCAPP")

    container "${params.SCAPP_image}"

    input:
    tuple val(sample), path(assemblyGraph), val(maxKmer), path(bam)

    output:
    tuple val("${sample}"), path("${sample}_plasmids.fasta.gz"), emit: plasmids, optional: true
    tuple file(".command.sh"), file(".command.out"), file(".command.err"), file(".command.log")
    tuple val("${sample}"), path("${sample}_plasmids_stats.tsv"), emit: plasmidsStats, optional: true

    shell:
    template("SCAPP.sh")
}


process pPlasClass {

    label 'medium'

    tag "$sample $binID"

    publishDir params.output, saveAs: { filename -> getOutput("${sample}", params.runid, "PlasClass", filename) }, \
        pattern: "{**.tsv,**.fasta.gz}"

    when params.steps.containsKey("plasmid") && params.steps.plasmid?.containsKey("PlasClass")

    container "${params.PlasClass_image}"

    input:
    tuple val(sample), val(binID), path(assembly)

    output:
    tuple val("${sample}"), val("${binID}"), path("${binID}_plasmids.fasta.gz"), emit: plasmids, optional: true
    tuple val("${sample}"), val("${binID}"), path("${binID}_probs.tsv"), emit: probabilities
    tuple val("${sample}_${binID}"), val("${output}"), val(params.LOG_LEVELS.INFO), file(".command.sh"), \
        file(".command.out"), file(".command.err"), file(".command.log"), emit: logs

    shell:
    output = getOutput("${sample}", params.runid, "PlasClass", "")
    template("PlasClass.sh")
}


process pPLSDB {

    label 'medium'

    tag "$sample $binID"

    publishDir params.output, saveAs: { filename -> getOutput("${sample}", params.runid, "PLSDB", filename) }, \
            pattern: "{**.tsv}"

    containerOptions " --user 1000:1000 --volume ${params.databases}:${params.databases} "

    when params.steps.containsKey("plasmid") && params.steps.plasmid.containsKey("PLSDB")

    container "${params.mash_image}"

    input:
    tuple val(sample), val(binID), path(plasmids)

    output:
    tuple val("${sample}"), val("${binID}"), path("${binID}.tsv"), emit: allHits
    tuple val("${sample}"), val("${binID}"), path("${binID}_kmerThreshold_*.tsv"), emit: filteredHitsMetadata
    tuple val("${sample}_${binID}"), val("${output}"), val(params.LOG_LEVELS.INFO), file(".command.sh"), \
        file(".command.out"), file(".command.err"), file(".command.log"), emit: logs

    shell:
    output = getOutput("${sample}", params.runid, "PLSDB", "")
    template("plsdb.sh")
}


workflow wPlasmidsList {
     take:
       samplesContigs
       samplesBins
     main:
       _wPlasmids(samplesContigs, samplesBins)
    emit:
      newPlasmids = _wPlasmids.out.newPlasmids
}


workflow _wPlasmids {
     take:
       samplesContigs
       samplesBins
     main:
       // search for new pladmids
       samplesBins | pSCAPP
       samplesContigs | pPlasClass

       SAMPLE_IDX = 0
       BIN_IDX = 1

       pSCAPP.out.plasmids | map { plasmids -> [plasmids[SAMPLE_IDX], "assembly", plasmids[BIN_IDX]] } \
	| set { newPlasmids }

       // check if there are known plasmids
       pPlasClass.out.plasmids \
	| mix(newPlasmids) | pPLSDB

       // store logs
       pPlasClass.out.logs | mix(pPLSDB.out.logs) | pDumpLogs
     emit:
       probabilities = pPlasClass.out.probabilities
       newPlasmids = newPlasmids

}
