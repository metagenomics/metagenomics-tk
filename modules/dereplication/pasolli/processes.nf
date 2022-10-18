
String getOutput(SAMPLE, RUNID, MODULE , TOOL, filename){

    def module = ""
    if(MODULE.isEmpty()){
       module = params.modules.binning.name + '/' +
          params.modules.binning.version.major + "." +
          params.modules.binning.version.minor + "." +
          params.modules.binning.version.patch;
    } else {
       module = MODULE
    }

    return SAMPLE + '/' + RUNID + '/' + module  +
          '/' + TOOL + '/' + filename
}

String getAggregatedOutput(RUNID, MODULE , TOOL, filename){

    def module = ""
    if(MODULE.isEmpty()){
       module = params.modules.binning.name + '/' +
          params.modules.binning.version.major + "." +
          params.modules.binning.version.minor + "." +
          params.modules.binning.version.patch;
    } else {
       module = MODULE
    }

    return "AGGREGATED" + '/' + RUNID + '/' + module  +
          '/' + TOOL + '/' + filename
}


process pMashSketchGenome {

    container "${params.mash_image}"

    label 'tiny'

    tag "Bin: ${binid}"

    when:
    run

    input:
    val(run)
    val(mashSketchParams)
    tuple val(binid), path("g.fa")

    output:
    path("${binid}.msh"), emit: sketch
    tuple env(GENOME_PATH), val("${binid}"), emit: stagedGenome
    tuple val("${binid}"), val("${binid}"), val(params.LOG_LEVELS.ALL), file(".command.sh"), \
	file(".command.out"), file(".command.err"), file(".command.log"), emit: logs

    shell:
    '''
    ln -s g.fa !{binid}
    mash sketch !{mashSketchParams} !{binid} -o !{binid}.msh
    GENOME_PATH=$(readlink -f g.fa)
    '''
}


process pMashPaste {

    container "${params.mash_image}"

    label 'large'

    publishDir params.output, mode: "${params.publishDirMode}", saveAs: { filename -> getAggregatedOutput(params.runid, "${module}", "${outputToolDir}", filename) }

    when:
    run

    input:
    val(run)
    tuple val(module), val(outputToolDir)
    path sketches, stageAs: 'sketch*.msh'

    output:
    path('final_sketch.msh'), emit: sketch
    tuple file(".command.sh"), file(".command.out"), file(".command.err"), file(".command.log"), emit: logs

    shell:
    '''
    mash paste final_sketch !{sketches}
    '''
}
