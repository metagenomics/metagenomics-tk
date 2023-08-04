nextflow.enable.dsl=2

include { getMetaflyeQualityParam } from '../assembly/ontAssembler'

def getOutput(SAMPLE, RUNID, TOOL, filename){
    return SAMPLE + '/' + RUNID + '/' + params.modules.assemblyHybrid.name + '/' +
          params.modules.assemblyHybrid.version.major + "." +
          params.modules.assemblyHybrid.version.minor + "." +
          params.modules.assemblyHybrid.version.patch +
          '/' + TOOL + '/' + filename
}

def timestamp = new java.util.Date().format( 'YYYYMMdd-HHmmss-SSS')

/*
*
* Metaflye processes long reads. Metaflye requires either --nano-raw, --nano-hq or --nano-hq and --read-error 0.03
* as input. The optimal parameter setting will be chosen by evaluating the median phred score.
*
*/
process pMetaflyeHybrid {

    label 'large'

    tag "Sample: $sample"

    publishDir params.output, mode: "${params.publishDirMode}", saveAs: { filename -> getOutput("${sample}", params.runid, "metaflye", filename) }

    when params?.steps.containsKey("assemblyHybrid") && params?.steps?.assemblyHybrid.containsKey("metaflye")

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
    
    shell:
    metaFlyeQualityParameter = params.steps?.assemblyHybrid?.metaflye.quality == "AUTO" ? \
	getMetaflyeQualityParam(Double.parseDouble(medianQuality)) : [quality: params.steps?.assemblyHybrid?.metaflye.quality, error: " "]
    quality = metaFlyeQualityParameter.quality
    error = metaFlyeQualityParameter.error
    '''
    ASSEMBLY_OUTPUT="!{sample}_contigs.fa.gz"
    HEADER_MAPPING_OUTPUT="!{sample}_contigs_header_mapping.tsv"

    flye !{quality} reads.fq.gz -o out --meta -t !{task.cpus} !{params.steps.assemblyHybrid.metaflye.additionalParams} !{error}

    # The following function modifies the assembly fasta headers according to the pattern: SAMPLEID_SEQUENCECOUNTER_SEQUENCEHASH
    transform.sh out/assembly.fasta ${ASSEMBLY_OUTPUT} ${HEADER_MAPPING_OUTPUT} !{sample} !{task.cpus}

    mv out/assembly_graph.gfa !{sample}_assembly_graph.gfa

    mv out/assembly_info.txt !{sample}_assembly_info.txt

    # get basic contig stats 
    paste -d$'\t' <(echo -e "SAMPLE\n!{sample}") <(seqkit stat -Ta ${ASSEMBLY_OUTPUT}) > !{sample}_contigs_stats.tsv
    '''
}


process pMetaspadesHybrid {

    label 'large'

    tag "$sample"

    publishDir params.output, mode: "${params.publishDirMode}", saveAs: { filename -> getOutput("${sample}", params.runid, "metaspades", filename) }

    when params?.steps.containsKey("assemblyHybrid") && params?.steps?.assemblyHybrid.containsKey("metaspades")

    container "${params.metaspades_image}"

    input:
    tuple val(sample), path(ontReads, stageAs: 'ontReads.fq.gz'), path(interleavedReads, stageAs: 'interleaved.fq.gz')

    output:
    tuple val("${sample}"), path("${sample}_contigs.fa.gz"), emit: contigs
    tuple val("${sample}"), path("${sample}_contigs_stats.tsv"), emit: contigsStats
    tuple val("${sample}"), path("${sample}_contigs.fastg"), env(maxKmer), emit: fastg, optional: true
    tuple file(".command.sh"), file(".command.out"), file(".command.err"), file(".command.log"), emit: logs

    shell:
    outputFastg = params?.steps?.assemblyHybrid?.metaspades?.fastg ? "TRUE" : "FALSE"
    template 'metaspades_hybrid.sh'
}

//TODO: check if generic method is available
process pBwa2Index {
    container "${params.bwa2_image}"
    label 'small'
    when params.steps.containsKey("assemblyHybrid") && params.steps.assemblyHybrid.containsKey("bwa2")
    input:
      tuple val(sample), path(contigs)
    output:
      tuple val("${sample}"), path('*.{0123,amb,ann,64,pac}'), emit: index
    shell:
      """
      bwa-mem2 index !{params.steps.assemblyHybrid.bwa2.additionalParams.bwa2_index} !{contigs}
      """
}

//TODO: check if generic method is available
process pMapBwa2 {
    label 'large'
    container "${params.samtools_bwa2_image}"
    when params.steps.containsKey("assemblyHybrid") && params.steps.assemblyHybrid.containsKey("bwa2")
    //publishDir params.output, mode: "${params.publishDirMode}", saveAs: { filename -> getOutput(params.runid ,"bwa2", filename) }
    input:
      tuple val(sample), path(reads), path(contigs), path(index, stageAs: "*") 

    output:
      tuple val("${sample}"), path("*bam"), path("*bai"), emit: alignment
      tuple val("${sample}"), val("${output}"), val(params.LOG_LEVELS.INFO), file(".command.sh"), \
        file(".command.out"), file(".command.err"), file(".command.log"), emit: logs

    shell:
      """
      bwa-mem2 mem !{params.steps.assemblyHybrid.bwa2.additionalParams.bwa2_mem} -p -t !{task.cpus} !{contigs} < !{reads} - \
        | samtools view -@ !{task.cpus} -S -b - \
        | samtools sort -l 9 -@ !{task.cpus} - > !{sample}.bam

      samtools index !{sample}.bam
      """
}

//TODO: use generic method from readmapping module (template needs fixes to be generic)
process pMinimap2IndexONT {

    container "${params.ubuntu_image}"

    label 'large'

    when params.steps.containsKey("assemblyHybrid") && params?.steps?.assemblyHybrid.containsKey("minimap")

    input:
      tuple val(sample), path(contigs)

    output:
      tuple val("${sample}"), path('seq.mmi'), emit: index
      

    shell:
      """
      minimap2 !{params.steps.assemblyHybrid.minimap.additionalParams.minimap_index} -x map-ont -d seq.mmi !{contigs}
      """
}

//TODO: use generic method from readmapping module (template needs fixes to be generic)
process pMapMinimap2ONT {
    label 'large'

    container "${params.samtools_bwa_image}"

    when params.steps.containsKey("assemblyHybrid") && params?.steps?.assemblyHybrid.containsKey("minimap")


    input:
      tuple val(sample), path(reads), path(index)

    output:
      tuple val("${sample}"), path("*sam.gz"), emit: alignment
      //tuple val("${sample}"), path("*bam"), path("*bam.bai"), emit: alignment
      tuple val("${sample}"), val("${output}"), val(params.LOG_LEVELS.INFO), file(".command.sh"), \
        file(".command.out"), file(".command.err"), file(".command.log"), emit: logs

    shell:
      """
      minimap2 !{params.steps.assemblyHybrid.minimap.additionalParams.minimap} -x map-ont -a !{index} -t !{task.cpus}  <(cat !{reads}) - \
            | gzip > !{sample}.sam.gz
      """
}

process pRacon {
    label 'large'

    container "${params.racon_image}"

    when params.steps.containsKey("assemblyHybrid") && params?.steps?.assemblyHybrid.containsKey("racon")

    publishDir params.output, mode: "${params.publishDirMode}", saveAs: { filename -> getOutput("${sample}", params.runid, "racon", filename) }

    input:
      tuple val(sample), path(reads), path(contigs), path(sam)

    output:
      tuple val("${sample}"), path("racon_polished.fasta"), emit: contigs
      tuple val("${sample}"), val("${output}"), val(params.LOG_LEVELS.INFO), file(".command.sh"), \
        file(".command.out"), file(".command.err"), file(".command.log"), emit: logs

    shell:
      """
      racon -m 8 -x -6 -g -8 -w 500 -t !{task.cpus} !{reads} !{sam} !{contigs} > racon_polished.fasta
      """
}

//TODO: medaka model is currently a parameter, might be better to add it to the input file in order to be able to process ont data of different flowcells
process pMedaka {
    label 'large'

    container "${params.medaka_image}"

    when params.steps.containsKey("assemblyHybrid") && params?.steps?.assemblyHybrid.containsKey("medaka")

    publishDir params.output, mode: "${params.publishDirMode}", saveAs: { filename -> getOutput("${sample}", params.runid, "medaka", filename) }

    input:
      tuple val(sample), path(reads), path(contigs)

    output:
      tuple val("${sample}"), path("medaka/consensus.fasta"), emit: contigs
      tuple val("${sample}"), val("${output}"), val(params.LOG_LEVELS.INFO), file(".command.sh"), \
        file(".command.out"), file(".command.err"), file(".command.log"), emit: logs

    shell:
      """
      medaka_consensus -i !{reads} -d !{contigs} -o medaka -t !{task.cpus} -m !{params.steps.assemblyHybrid.medaka.model}
      """
}

process pPilon {
    label 'small'

    container "${params.pilon_image}"

    when params.steps.containsKey("assemblyHybrid") && params?.steps?.assemblyHybrid.containsKey("pilon")

    publishDir params.output, mode: "${params.publishDirMode}", saveAs: { filename -> getOutput("${sample}", params.runid, "pilon", filename) }

    input:
      tuple val(sample), path(reads), path(contigs), path(bam), path(bai)

    output:
      tuple val("${sample}"), path("pilon.fasta"), emit: contigs
      tuple val("${sample}"), val("${output}"), val(params.LOG_LEVELS.INFO), file(".command.sh"), \
        file(".command.out"), file(".command.err"), file(".command.log"), emit: logs

    shell:
      """
      pilon --genome !{contigs} --fix all --changes --frags !{bam} --threads !{task.cpus} --output pilon
      """
}


/*
 * Input: 
 *  - List of the format [SAMPLE, ONT Reads, CONTIGS]
 *
 * Output is of the format [SAMPLE, CONTIGS].
 * 
 * Runs minimap2 for mapping of ONT reads to assembly and runs racon for polishing
 */
workflow _wRaconPolishing {
    take:
       polishingInput
    main:
       SAMPLE_IDX = 0
       ONT_IDX = 1
       CONTIGS_IDX = 2

       //create index
       polishingInput | map {it -> [it[SAMPLE_IDX], it[CONTIGS_IDX]] } | pMinimap2IndexONT
       //run minimap mapping
       polishingInput | map {it -> [it[SAMPLE_IDX], it[ONT_IDX]] } | join(pMinimap2IndexONT.out.index) | pMapMinimap2ONT
       //run racon
       polishingInput | join(pMapMinimap2ONT.out.alignment) | pRacon
    emit:
      contigs = pRacon.out.contigs
}


/*
 * Input: 
 *  - List of the format [SAMPLE, Illumina Reads, CONTIGS]
 *
 * Output is of the format [SAMPLE, CONTIGS].
 * 
 * Runs bwa2 for mapping of ONT reads to assembly and runs pilon for polishing
 */
workflow _wPilonPolishing {
    take:
       polishingInput
    main:
       SAMPLE_IDX = 0
       PAIRED_END_IDX = 1
       CONTIGS_IDX = 2

       //create index
       polishingInput | map {it -> [it[SAMPLE_IDX], it[CONTIGS_IDX]] } | pBwa2Index
       //run minimap mapping
       polishingInput | map {it -> [it[SAMPLE_IDX], it[PAIRED_END_IDX], it[CONTIGS_IDX]] } | join(pBwa2Index.out.index) | pMapBwa2
       //run Pilon
       polishingInput | join(pMapBwa2.out.alignment) | pPilon
    emit:
      contigs = pPilon.out.contigs
}

/*
 * Input: 
 *  - List of the format [SAMPLE, ONT Reads, Reads Paired, Reads Unpaired] 
 *
 * Output is of the format [SAMPLE, CONTIGS] for contigs and [SAMPLE, fastg] for fastq files.
 * 
 */
workflow wHybridAssemblyList {
    take:
       readsList
    main:
       _wHybridAssembly(readsList)
    emit:
      contigs = _wHybridAssembly.out.contigs
      graph = _wHybridAssembly.out.graph
      contigsStats = _wHybridAssembly.out.graph
}

/*
*
* Input:
*  - List of the format [SAMPLE, ONT Reads, Reads Paired, Reads Unpaired]
*/
workflow _wHybridAssembly {
     take:
       readsList
     main:
       SAMPLE_IDX = 0
       ONT_IDX = 1
       PAIRED_END_IDX = 2
       UNPAIRED_END_IDX = 3
       ONT_MEDIAN_QUALITY_IDX = 4
 
       //start metaspades hybrid assembly
       readsList | map { seq -> [seq[SAMPLE_IDX], seq[ONT_IDX], seq[PAIRED_END_IDX]] } | pMetaspadesHybrid
       
       //TODO: check if this is needed:
       if(params.summary){
         pMetaspadesHybrid.out.contigsStats | collectFile(newLine: false, keepHeader: true, storeDir: params.output + "/summary/" ){ item ->
           [ "contigs_stats.tsv", item[1].text ]
         }
        }


       //start metaflye ONT assembly:
       readsList | map { seq -> [seq[SAMPLE_IDX], seq[ONT_IDX], seq[ONT_MEDIAN_QUALITY_IDX]] } | pMetaflyeHybrid
       //run Racon polishing:
       readsList | map { seq -> [seq[SAMPLE_IDX], seq[ONT_IDX]] } | join(pMetaflyeHybrid.out.contigs) | _wRaconPolishing
       //run medaka polishing:
       readsList | map { seq -> [seq[SAMPLE_IDX], seq[ONT_IDX]] } | join(_wRaconPolishing.out.contigs) |  pMedaka
       //run pilon polishing (with bwa2 and minimap2 mappings)
       readsList | map { seq -> [seq[SAMPLE_IDX], seq[PAIRED_END_IDX]] } | join(pMedaka.out.contigs) | _wPilonPolishing

       //TODO: additional outputs for binning might be necessary?
       pMetaspadesHybrid.out.contigs | mix(_wPilonPolishing.out.contigs) | set { contigs }
       
       //TODO: fastg output and gfa mixed together - problem?
       pMetaspadesHybrid.out.fastg | mix(pMetaflyeHybrid.out.graph) | set { graph }
       
       pMetaspadesHybrid.out.contigsStats | mix(pMetaflyeHybrid.out.contigsStats) | set { contigsStats }

    emit:
      contigs = contigs
      graph = graph 
      contigsStats = contigsStats
}
