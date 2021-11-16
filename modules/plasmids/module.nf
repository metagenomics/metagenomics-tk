nextflow.enable.dsl=2


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
    tuple val("${sample}"), path("${sample}_plasmids.fasta.gz"), emit: plasmids
    tuple file(".command.sh"), file(".command.out"), file(".command.err"), file(".command.log")
    tuple val("${sample}"), path("${sample}_plasmids_stats.tsv"), emit: plasmidsStats

    shell:
    template("SCAPP.sh")
}


process pPlasClass {

    label 'medium'

    tag "$binID"

    publishDir params.output, saveAs: { filename -> getOutput("${sample}", params.runid, "PlasClass", filename) }

    when params.steps.containsKey("plasmid") && params.steps.plasmid?.containsKey("PlasClass")

    container "${params.PlasClass_image}"

    input:
    tuple val(sample), val(binID), path(assembly)

    output:
    tuple val("${sample}"), val("${binID}"), path("${binID}_plasmids.fasta.gz"), emit: plasmids, optional: true
    tuple val("${sample}"), val("${binID}"), path("${binID}_probs.tsv"), emit: probabilities
    tuple file(".command.sh"), val("${binID}"), file(".command.out"), file(".command.err"), file(".command.log")

    shell:
    template("PlasClass.sh")
}


process pPLSDB {

    label 'small'

    tag "$binID"

    publishDir params.output, saveAs: { filename -> getOutput("${sample}", params.runid, "PLSDB", filename) }

    when params.steps.containsKey("plasmid") && params.steps.plasmid?.containsKey("PLSDB")

    container "${params.mash_image}"

    input:
    tuple val(sample), val(binID), path(plasmids)

    output:
    tuple file(".command.sh"), file(".command.out"), file(".command.err"), file(".command.log")
    tuple val("${sample}"), val("${binID}"), path("${binID}.tsv"), emit: allHits
    tuple val("${sample}"), val("${binID}"), path("${binID}_kmerThreshold_*.tsv"), emit: filteredHitsMetadata

    shell:
    template("plsdb.sh")
}


workflow wPlasmidsList {
     take:
       samplesAssembly
       samplesBins
     main:
       _wPlasmids(samplesAssembly, samplesBins)
    emit:
      newPlasmids = _wPlasmids.out.newPlasmids
}


workflow _wPlasmids {
     take:
       samplesAssembly
       samplesBins
     main:
       samplesBins | pSCAPP
       samplesAssembly | pPlasClass
       pPlasClass.out.plasmids | mix(pSCAPP.out.plasmids | map { plasmids -> [plasmids[0], "assembly", plasmids[1]] }) | pPLSDB
       pSCAPP.out.plasmids | map { plasmids -> [plasmids[0], "plasmid", plasmids[1]] } \
	| set { newPlasmids }
     emit:
       probabilities = pPlasClass.out.probabilities
       newPlasmids = newPlasmids
}
