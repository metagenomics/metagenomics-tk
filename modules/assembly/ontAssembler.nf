include { wSaveSettingsList } from '../config/module'


def getOutput(SAMPLE, RUNID, TOOL, filename){
    return SAMPLE + '/' + RUNID + '/' + params.modules.assemblyONT.name + '/' +
          params.modules.assemblyONT.version.major + "." +
          params.modules.assemblyONT.version.minor + "." +
          params.modules.assemblyONT.version.patch +
          '/' + TOOL + '/' + filename
}


def getMetaflyeQualityParam(medianQuality) {
   if(medianQuality > 20){
     return [ quality: " --nano-hq ", error: " --read-error 0.03 "]
   }
   if(medianQuality > 13){
     return [ quality: " --nano-hq ", error: ""]
   } else {
     return [ quality: " --nano-raw ", error: ""]
   }
}

/*
*
* Metaflye processes long reads. Metaflye requires either --nano-raw, --nano-hq or --nano-hq and --read-error 0.03
* as input. The optimal parameter setting will be chosen by evaluating the median phred score.
*
*/
process pMetaflye {

    label 'highmemLarge'

    tag "Sample: $sample"

    publishDir params.output, mode: "${params.publishDirMode}", saveAs: { filename -> getOutput("${sample}", params.runid, "metaflye", filename) }

    when params?.steps.containsKey("assemblyONT") && params?.steps?.assemblyONT.containsKey("metaflye")

    container "${params.metaflye_image}"

    input:
    tuple val(sample), path(reads, stageAs: 'reads.fq.gz'), val(medianQuality)

    output:
    tuple val("${sample}"), path("${sample}_contigs.fa.gz"), emit: contigs
    tuple val("${sample}"), path("${sample}_assembly_info.txt"), emit: info
    tuple val("${sample}"), path("${sample}_contigs_header_mapping.tsv"), emit: headerMapping
    tuple val("${sample}"), path("${sample}_assembly_graph.gfa"), emit: graph
    tuple val("${sample}"), path("${sample}_contigs_stats.tsv"), emit: contigsStats
    tuple file(".command.sh"), file(".command.out"), file(".command.err"), file(".command.log")
    
    script:
    metaFlyeQualityParameter = params.steps?.assemblyONT?.metaflye.quality == "AUTO" ? \
	getMetaflyeQualityParam(Double.parseDouble(medianQuality)) : [quality: params.steps?.assemblyONT?.metaflye.quality, error: " "]
    quality = metaFlyeQualityParameter.quality
    error = metaFlyeQualityParameter.error
    """
    ASSEMBLY_OUTPUT="${sample}_contigs.fa.gz"
    HEADER_MAPPING_OUTPUT="${sample}_contigs_header_mapping.tsv"

    flye ${quality} reads.fq.gz -o out --meta -t ${task.cpus} ${params.steps.assemblyONT.metaflye.additionalParams} ${error}

    # The following function modifies the assembly fasta headers according to the pattern: SAMPLEID_SEQUENCECOUNTER_SEQUENCEHASH
    transform.sh out/assembly.fasta \${ASSEMBLY_OUTPUT} \${HEADER_MAPPING_OUTPUT} ${sample} ${task.cpus}

    mv out/assembly_graph.gfa ${sample}_assembly_graph.gfa

    mv out/assembly_info.txt ${sample}_assembly_info.txt

    # get basic contig stats 
    paste -d\$'\t' <(echo -e "SAMPLE\n${sample}") <(seqkit stat -Ta \${ASSEMBLY_OUTPUT}) > ${sample}_contigs_stats.tsv
    """
}

/*
*
* metaMDBG only runs for flowcells > 10.4.1 
*
*/
process pMetaMDBG {
    label 'highmemLarge'
    tag "Sample: $sample"

    container "quay.io/biocontainers/metamdbg:1.4--h3be2455_0"

    publishDir params.output, mode: "${params.publishDirMode}", saveAs: { filename -> getOutput("${sample}", params.runid, "metaMDBG", filename) }

    when:
    params?.steps?.containsKey("assemblyONT") && params?.steps?.assemblyONT?.containsKey("metaMDBG")

    input:
    tuple val(sample), path(reads, stageAs: 'reads.fq.gz'), val(medianQuality)

    output:
    tuple val("${sample}"), path("${sample}_contigs.fa.gz"), emit: contigs
    tuple val("${sample}"), path("${sample}_contigs_header_mapping.tsv"), emit: headerMapping
    tuple val("${sample}"), path("${sample}_assembly_graph.gfa"), emit: graph
    tuple val("${sample}"), path("${sample}_contigs_stats.tsv"), emit: contigsStats
    tuple file(".command.sh"), file(".command.out"), file(".command.err"), file(".command.log")

    shell:
    // Erlaube zusätzliche Parameter und einen konfigurierbaren k-Wert für den Graphen (Standard: 21)
    additionalParams = params.steps.assemblyONT.metaMDBG.additionalParams.metaMDBG ?: ""
    kVal = params.steps.assemblyONT.metaMDBG.additionalParams.kVal ?: "21"
    
    '''
    ASSEMBLY_OUTPUT="!{sample}_contigs.fa.gz"
    HEADER_MAPPING_OUTPUT="!{sample}_contigs_header_mapping.tsv"

    # 1. Assembly
    metaMDBG asm --out-dir out --in-ont reads.fq.gz --threads !{task.cpus} !{additionalParams}
    
    # 2. If k_val is set use that for graph generation, otherwise take the largest
    K_VAL="!{kVal}"

    if [ -z "$K_VAL" ]; then
    
        #find biggest k:
        MAX_K=$(metaMDBG gfa --assembly-dir out --k 0 2>&1 | grep '-' | awk -F'-' '{print $2}' | awk '{print $1}' | sort -n | tail -n 1)    
        
        if [ -n "$MAX_K" ]; then
            K_VAL=$MAX_K
        fi
    fi

    # 3. Graphen mit ermitteltem K_VAL generieren
    metaMDBG gfa --assembly-dir out --k $K_VAL --threads !{task.cpus}

    # Uncompress contigs for transform.sh
    gunzip out/contigs.fasta.gz

    # Header Transformation
    transform.sh out/contigs.fasta ${ASSEMBLY_OUTPUT} ${HEADER_MAPPING_OUTPUT} !{sample} !{task.cpus}

    # Find and rename assembly graph
    GFA_FILE=$(ls out/*.gfa 2>/dev/null | head -n 1)
    if [ -n "$GFA_FILE" ]; then
        mv "$GFA_FILE" !{sample}_assembly_graph.gfa
    else
        # Fallback, falls die Graph-Generierung fehlschlägt oder übersprungen wird
        touch !{sample}_assembly_graph.gfa
    fi

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
      graph = _wOntAssembly.out.graph
      mapping = _wOntAssembly.out.mapping
      headerMapping = _wOntAssembly.out.mapping
      info = _wOntAssembly.out.info
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
             | map { it -> [ it.SAMPLE, it.READS, file("NOT_SET")]} | set { reads }
       
       _wOntAssembly(reads)

       SAMPLE_IDX = 0
       wSaveSettingsList(reads | map { it -> it[SAMPLE_IDX] })
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
       readsList | pMetaMDBG 
       pMetaflye.out.contigs | mix(pMetaMDBG.out.contigs) | set { contigs }
       pMetaflye.out.graph | mix(pMetaMDBG.out.graph) | set { graph }
       pMetaflye.out.headerMapping | mix(pMetaMDBG.out.headerMapping) | set { headerMapping }

     emit:
       contigs = contigs
       graph = graph
       mapping = headerMapping
       info = pMetaflye.out.info
}
