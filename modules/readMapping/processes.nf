


process pMinimap2Index {
    container "${params.ubuntu_image}"
    label 'large'

    when:
    run

    input:
      val(run)
      val(seqType)
      path(representatives)
    output:
      path('seq.mmi')
    shell:
      """
      minimap2 !{params.steps.readMapping.minimap.additionalParams.minimap_index} -x !{seqType} -d seq.mmi !{representatives}
      """
}

process pMapMinimap2 {
    label 'large'
    container "${params.samtools_bwa_image}"
    when:
    run
    publishDir params.output, mode: "${params.publishDirMode}", saveAs: { filename -> getOutput("${sampleID}", params.runid ,"minimap", filename) }
    input:
      val(run)
      tuple val(sampleID), path(sample, stageAs: "sample*"), path(index)

    output:
      tuple val("${sampleID}"), path("*bam"), path("*bam.bai"), emit: alignment
      tuple file(".command.sh"), file(".command.out"), file(".command.err"), file(".command.log"), emit: logs
    shell:
    template('minimap2.sh')
}


