nextflow.enable.dsl=2

MODULE="assembly"
VERSION="1.0.0"
def getOutput(SAMPLE, RUNID, TOOL, filename){
    return SAMPLE + '/' + RUNID + '/' + MODULE + '/' + VERSION + '/' + TOOL + '/' + filename
}


process pMegahit {

    label 'large'

    tag "$sample"

    publishDir params.output, saveAs: { filename -> getOutput("${sample}", params.runid, "megahit", filename) }

    errorStrategy 'ignore'

    when params?.steps.containsKey("assembly") && params?.steps?.assembly.containsKey("megahit")

    container "vout/megahit:${params.megahit_tag}"

    input:
    tuple val(sample), path(fastqs, stageAs: 'reads.fq.gz')

    output:
    tuple val("${sample}"), path("${sample}_contigs.fa.gz"), emit: contigs
    tuple val("${sample}"), path("${sample}_contigs_stats.tsv"), emit: contigs_stats

    shell:
    template 'megahit.sh'
}



/*
 * Takes a list as input with the format [SAMPLE, READS]
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
     take:
       readsTable
     main:
       readsTable | splitCsv(sep: '\t', header: true) \
             | map { it -> [ it.SAMPLE, it.READS]} \
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
         pMegahit.out.contigs_stats | collectFile(newLine: false, keepHeader: true, storeDir: params.output + "/summary/" ){ item ->
           [ "contigs_stats.tsv", item[1].text ]
         }
       }

    emit:
      contigs = pMegahit.out.contigs
}
