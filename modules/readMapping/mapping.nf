nextflow.enable.dsl=2


def getOutput(SAMPLE, RUNID, TOOL, filename){
    return SAMPLE + '/' + RUNID + '/' + params.modules.readMapping.name + '/' + 
         params.modules.readMapping.version.major + "." + 
         params.modules.readMapping.version.minor + "." + 
         params.modules.readMapping.version.patch +
         '/' + TOOL + '/' + filename
}

process pBwaIndex {
    container "${params.bwa_image}"
    label 'large'
    when params.steps.containsKey("readMapping") && params.steps.readMapping.containsKey("bwa")
    input:
      path(representatives)
    output:
      path('*.{amb,ann,bwt,pac,sa,fa}')
    shell:
      """
      bwa index !{params.steps.readMapping.bwa.additionalParams.bwa_index} !{representatives}
      """
}


process pMinimap2Index {
    container "${params.ubuntu_image}"
    label 'large'
    when params.steps.containsKey("readMapping") && params.steps.readMapping.containsKey("minimap2")
    input:
      path(representatives)
    output:
      path('ont.mmi')
    shell:
      """
      minimap2 !{params.steps.readMapping.minimap2.additionalParams.minimap2_index} -x map-ont -d ont.mmi !{representatives}
      """
}


process pMapMinimap2 {
    label 'large'
    container "${params.samtools_bwa_image}"
    when params.steps.containsKey("readMapping") && params.steps.readMapping.containsKey("minimap2")
    publishDir params.output, mode: "${params.publishDirMode}", saveAs: { filename -> getOutput("${sampleID}", params.runid ,"minimap2", filename) }
    input:
      tuple val(sampleID), path(sample), val(mode), path(representatives_fasta), path(x, stageAs: "*") 
    output:
      tuple val("${sampleID}"), path("*bam"), path("*bam.bai"), emit: alignment
      tuple file(".command.sh"), file(".command.out"), file(".command.err"), file(".command.log"), emit: logs
    shell:
    template('minimap2.sh')
}


process pMapBwa {
    label 'large'
    container "${params.samtools_bwa_image}"
    when params.steps.containsKey("readMapping")
    publishDir params.output, mode: "${params.publishDirMode}", saveAs: { filename -> getOutput("${sampleID}", params.runid ,"bwa", filename) }
    input:
      tuple val(sampleID), path(sample), val(mode), path(representatives_fasta), path(x, stageAs: "*") 
    output:
      tuple val("${sampleID}"), path("*bam"), emit: alignment
      tuple file(".command.sh"), file(".command.out"), file(".command.err"), file(".command.log"), emit: logs
    shell:
    MODE = mode == "paired" ? " -p " : ""
    template('bwa.sh')
}


process pMergeAlignment {
    when params.steps.containsKey("readMapping")
    label 'tiny'
    container 'quay.io/biocontainers/samtools:1.12--h9aed4be_1'
    publishDir params.output, mode: "${params.publishDirMode}", saveAs: { filename -> getOutput("${sample}", params.runid, "coverm", filename) }
    input:
      tuple val(sample), file("alignment?.bam")
    output:
      tuple val("${sample}"), path("${sample}.bam"), path("*bam.bai"), emit: alignmentIndex
      tuple file(".command.sh"), file(".command.out"), file(".command.err"), file(".command.log")
    shell:
    """
    samtools merge !{sample}.bam alignment*.bam
    samtools index !{sample}.bam
    """
}


process pCovermCount {
    when params.steps.containsKey("readMapping")
    label 'small'
    publishDir params.output, mode: "${params.publishDirMode}", saveAs: { filename -> getOutput("${sample}", params.runid, "coverm", filename) }
    input:
      tuple val(sample), file(mapping), file(index), file(list_of_representatives), val(medianQuality)
    output:
      tuple val("${sample}"), path("${sample}_out/mean.tsv"), emit: mean
      tuple val("${sample}"), path("${sample}_out/trimmed_mean.tsv"), emit: trimmedMean
      tuple val("${sample}"), path("${sample}_out/count.tsv"), emit: count
      tuple val("${sample}"), path("${sample}_out/rpkm.tsv"), emit: rpkm
      tuple val("${sample}"), path("${sample}_out/tpm.tsv"), emit: tpm
      tuple file(".command.sh"), file(".command.out"), file(".command.err"), file(".command.log")
    shell:
    DO_NOT_ESTIMATE_QUALITY = -1 
    MEDIAN_QUALITY=Double.parseDouble(medianQuality)
    percentIdentity = MEDIAN_QUALITY != DO_NOT_ESTIMATE_QUALITY ? \
	" --min-read-percent-identity "+Utils.getMappingIdentityParam(MEDIAN_QUALITY) : " "
    template('coverm.sh')
}


/*
*
* wListReadMappingBwa maps a set of fastq samples against a set of genomes.
* Input Parameters
*
*  * Columns of the mags file must be SAMPLE for sample identifier and READS for the path to the fastq files.
*
*/
workflow wFileReadMappingBwa {
   main:
     GENOMES_PATH_IDX = 0
     Channel.from(file(params?.steps?.readMapping?.mags)) \
       | splitCsv(sep: '\t', header: false, skip: 1) \
       | map { it -> file(it[GENOMES_PATH_IDX]) } \
       | set {genomesList}

     Channel.from(file(params?.steps?.readMapping?.samples)) \
       | splitCsv(sep: '\t', header: true)\
       | map { it -> [it.SAMPLE, it.READS] } | set {samples}

     _wReadMappingBwa(samples, Channel.empty(), genomesList)
   emit:
     trimmedMean = _wReadMappingBwa.out.trimmedMean
}

/*
*
* wListReadMappingBwa maps a set of fastq samples against a set of genomes.
* Input Parameters
*
*  * Entries in the samples channel consist of Maps of the format [SAMPLE:test1, READS:/path/to/interleaved/fastq/test1_interleaved.qc.fq.gz]. 
*
*  * Entries of the genomes channel are paths to the genomes. 
*
*/
workflow wListReadMappingBwa {
   take:
     samplesONT
     ontMedianQuality
     samplesPaired
     samplesSingle
     genomes
   main:
     _wReadMappingBwa(samplesONT, ontMedianQuality, samplesPaired, samplesSingle, genomes)
   emit:
     trimmedMean = _wReadMappingBwa.out.trimmedMean
}


workflow _wReadMappingBwa {
   take:
     samplesONT
     ontMedianQuality
     samplesPaired
     samplesSingle
     genomes
   main:
     // Create temporary directory for merged fasta files
     genomesTempDir = params.tempdir + "/genomesToBeMerged"
     file(genomesTempDir).mkdirs()

     // Concatenate genomes 
     genomes | collectFile(tempDir: genomesTempDir){ item -> [ "mergedGenomes.fasta", item.text ] } \
      | set { genomesMerged }

     // Create BWA and Minimap index of all genomes
     BWA_INDEX_IDX=0
     GENOMES_IDX=1
     genomesMerged | pBwaIndex | set {index} 

     genomesMerged | pMinimap2Index | set { ontIndex }

     // combine index with every sample
     index | map{ bwaIndex -> [bwaIndex]} \
      | combine(genomesMerged) \
      | map{ it -> ["paired", it[GENOMES_IDX], it[BWA_INDEX_IDX]] } \
      | set {pairedIndex}
     samplesPaired | combine(pairedIndex) | set {paired}

     index | map{ bwaIndex -> [bwaIndex]} \
      | combine(genomesMerged) \
      | map{ it -> ["single", it[GENOMES_IDX], it[BWA_INDEX_IDX]] } \
      | set {singleIndex}
     samplesSingle | combine(singleIndex) | set {single}

     ontIndex | map{ bwaIndex -> [bwaIndex]} \
      | combine(genomesMerged)  \
      | map{ it -> ["ONT", it[GENOMES_IDX], it[BWA_INDEX_IDX]] } \
      | set {ontLongReadIndex}
     samplesONT | combine(ontLongReadIndex) | set {ont}

     // Paired and single reads should be mapped back
     pMapBwa(paired | mix(single))

     // Map ONT data
     pMapMinimap2(ont)
 
     // The resulting alignments (bam files) should merged if single and paired read alignments exist
     SAMPLE_NAME_IDX=0
     pMapBwa.out.alignment | groupTuple(by: SAMPLE_NAME_IDX) | pMergeAlignment

     DO_NOT_ESTIMATE_IDENTITY = "-1"
     pMergeAlignment.out.alignmentIndex | combine(genomes | map {it -> file(it)} \
      | toList() | map { it -> [it]})  \
      | combine(Channel.value(DO_NOT_ESTIMATE_IDENTITY)) \
      | set { covermBWAInput }

     pMapMinimap2.out.alignment | combine(genomes | map {it -> file(it)} \
      | toList() | map { it -> [it]}) \
      | join(ontMedianQuality, by: SAMPLE_NAME_IDX) | set { covermMinimapInput }

     covermBWAInput | mix(covermMinimapInput) | pCovermCount
   emit:
     trimmedMean = pCovermCount.out.trimmedMean
}
