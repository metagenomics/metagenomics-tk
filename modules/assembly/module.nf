nextflow.enable.dsl=2


def getOutput(SAMPLE, RUNID, TOOL, filename){
    return SAMPLE + '/' + RUNID + '/' + params.modules.assembly.name + '/' +
          params.modules.assembly.version.major + "." +
          params.modules.assembly.version.minor + "." +
          params.modules.assembly.version.patch +
          '/' + TOOL + '/' + filename
}


process pMegahit {

    label 'large'

    tag "$sample"

    publishDir params.output, saveAs: { filename -> getOutput("${sample}", params.runid, "megahit", filename) }

    when params?.steps.containsKey("assembly") && params?.steps?.assembly.containsKey("megahit")

    container "${params.megahit_image}"

    input:
    tuple val(sample), path(interleavedReads, stageAs: 'interleaved.fq.gz'), path(unpairedReads)

    output:
    tuple val("${sample}"), path("${sample}_contigs.fa.gz"), emit: contigs
    tuple val("${sample}"), path("${sample}_contigs_stats.tsv"), emit: contigsStats
    tuple val("${sample}"), path("${sample}_contigs.fastg"), env(maxKmer), optional: true, emit: fastg
    tuple file(".command.sh"), file(".command.out"), file(".command.err"), file(".command.log")

    shell:
    includeUnpairedReads = unpairedReads.name != "NOT_SET" ? " -r ${unpairedReads} " : ''
    convertToFastg = params.steps.assembly.megahit.fastg ? "TRUE" : "FALSE"
    template 'megahit.sh'
}



/*
 * Takes a list as input with the format [SAMPLE, READS_PAIRED, READS_UNPAIRED]
 * Output is of the format [SAMPLE, CONTIGS]
 * 
 */
workflow wAssemblyList {
     take:
       readsList
     main:
       readsList | _wAssembly
    emit:
      contigs = _wAssembly.out.contigs
      fastg = _wAssembly.out.fastg
}



/*
 * Takes a tab separated file of files containing reads as input and produces assembly results.
 * Input file with columns seperated by tabs:
 * SAMPLE and READS
 * 
 * Output is of the format [SAMPLE, CONTIGS]
 * 
 */
workflow wAssemblyFile {
    main:
       Channel.from(file(params.steps.assembly.input)) | splitCsv(sep: '\t', header: true) \
             | map { it -> [ it.SAMPLE, it.READS, file("NOT_SET")]} \
             | _wAssembly
    emit:
      contigs = _wAssembly.out.contigs
}



/*
*
* Input: List of the format [SAMPLE, READS] 
*
*/
workflow _wAssembly {
     take:
       readsList
     main:
       readsList | pMegahit
       if(params.summary){
         pMegahit.out.contigsStats | collectFile(newLine: false, keepHeader: true, storeDir: params.output + "/summary/" ){ item ->
           [ "contigs_stats.tsv", item[1].text ]
         }
       }
    emit:
      contigs = pMegahit.out.contigs
      fastg = pMegahit.out.fastg
}
