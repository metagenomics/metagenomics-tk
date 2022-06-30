nextflow.enable.dsl=2


def getOutput(SAMPLE, RUNID, TOOL, filename){
    return SAMPLE + '/' + RUNID + '/' + params.modules.assemblyONT.name + '/' +
          params.modules.assemblyONT.version.major + "." +
          params.modules.assemblyONT.version.minor + "." +
          params.modules.assemblyONT.version.patch +
          '/' + TOOL + '/' + filename
}


process pMetaflye {

    label 'large'

    tag "Sample: $sample"

    publishDir params.output, saveAs: { filename -> getOutput("${sample}", params.runid, "metaflye", filename) }

    when params?.steps.containsKey("assemblyONT") && params?.steps?.assemblyONT.containsKey("metaflye")

    container "${params.metaflye_image}"

    input:
    tuple val(sample), path(reads, stageAs: 'reads.fq.gz')

    output:
    tuple val("${sample}"), path("${sample}_contigs.fa.gz"), emit: contigs
    tuple val("${sample}"), path("${sample}_contigs_stats.tsv"), emit: contigsStats
    tuple file(".command.sh"), file(".command.out"), file(".command.err"), file(".command.log")

    shell:
    '''
    ASSEMBLY_OUTPUT="!{sample}_contigs.fa.gz"
    flye --nano-raw reads.fq.gz -o out --meta -t !{task.cpus} !{params.steps.assemblyONT.metaflye.additionalParams}

    # The following function modifies the assembly fasta headers according to the pattern: SAMPLEID_SEQUENCECOUNTER_SEQUENCEHASH
    transform.sh out/assembly.fasta ${ASSEMBLY_OUTPUT} !{sample} !{task.cpus}

    # get basic contig stats 
    paste -d$'\t' <(echo -e "SAMPLE\n!{sample}") <(seqkit stat -Ta ${ASSEMBLY_OUTPUT}) > !{sample}_contigs_stats.tsv
    '''
}


/*
 * Takes a list as input with the format [SAMPLE, READS]
 * Output is of the format [SAMPLE, CONTIGS]
 * 
 */
workflow wOntAssemblyList {
     take:
       readsList
     main:
       readsList | _wOntAssembly
    emit:
      contigs = _wOntAssembly.out.contigs
}



/*
 * Takes a tab separated file of files containing reads as input and produces assembly results.
 * Input file with columns seperated by tabs:
 * SAMPLE and READS
 * 
 * Output is of the format [SAMPLE, CONTIGS]
 * 
 */
workflow wOntAssemblyFile {
    main:
       Channel.from(file(params.steps.assemblyONT.input)) | splitCsv(sep: '\t', header: true) \
             | map { it -> [ it.SAMPLE, it.READS, file("NOT_SET")]} \
             | _wOntAssembly
    emit:
      contigs = _wOntAssembly.out.contigs
}



/*
*
* Input: List of the format [SAMPLE, READS] 
*
*/
workflow _wOntAssembly {
     take:
       readsList
     main:
       readsList | pMetaflye 

       if(params.summary){
         pMetaflye.out.contigsStats | mix(pMetaflye.out.contigsStats) \
	  | collectFile(newLine: false, keepHeader: true, storeDir: params.output + "/summary/" ){ item ->
           [ "contigs_stats.tsv", item[1].text ]
          }
       }
            
     emit:
       contigs = pMetaflye.out.contigs
}
