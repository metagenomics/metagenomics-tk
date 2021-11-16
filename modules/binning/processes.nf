

def getOutput(SAMPLE, RUNID, TOOL, filename){
    return SAMPLE + '/' + RUNID + '/' + params.modules.binning.name + '/' +
          params.modules.binning.version.major + "." +
          params.modules.binning.version.minor + "." +
          params.modules.binning.version.patch +
          '/' + TOOL + '/' + filename
}

process pGetBinStatistics {

    container "${params.samtools_image}"

    tag "$sample"

    publishDir params.output, saveAs: { filename -> getOutput("${sample}", params.runid, "${binner}", filename) }

    label 'tiny'

    input:
    tuple val(sample), path(binContigMapping), path(bam), val(binner), path(bins)

    output:
    tuple val("${sample}"), file("${sample}_contigs_depth.tsv"), optional: true, emit: contigsDepth
    tuple val("${sample}"), file("${sample}_bins_stats.tsv"), optional: true, emit: binsStats
    tuple file(".command.sh"), file(".command.out"), file(".command.err"), file(".command.log")

    shell:
    template 'binStats.sh'
}
