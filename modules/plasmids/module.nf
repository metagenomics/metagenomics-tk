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


workflow wPlasmidsList {
     take:
       samples
     main:
       samples | _wPlasmids
    emit:
      plasmids = _wPlasmids.out.plasmids
}


workflow _wPlasmids {
     take:
       samples
     main:
       samples | pSCAPP 
     emit:
       plasmids = pSCAPP.out.plasmids
}
