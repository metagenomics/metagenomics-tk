
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

    fair true

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

    label 'highmemLarge'

    publishDir params.output, mode: "${params.publishDirMode}", \
	saveAs: { filename -> getAggregatedOutput(params.runid, "${module}", "${outputToolDir}", filename) }, \
        pattern: "{**.msh}"

    when:
    run

    input:
    val(run)
    tuple val(module), val(outputToolDir)
    path sketches, stageAs: 'sketch*.msh'

    output:
    path('final_sketch_*.msh'), emit: sketch
    tuple env(FILE_ID), val("${output}"), val(params.LOG_LEVELS.ALL), file(".command.sh"), \
	file(".command.out"), file(".command.err"), file(".command.log"), emit: logs

    shell:
    output = getAggregatedOutput(params.runid, module,"${outputToolDir}" , "")
    '''
    FILE_ID=$(mktemp XXXXXXXX)
    mash paste final_sketch_${FILE_ID} !{sketches}
    '''
}
