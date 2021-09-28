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
    when params.steps.containsKey("readMapping")
    input:
      path(representatives)
    output:
      path('*.{amb,ann,bwt,pac,sa,fa}')
    shell:
      """
      bwa index !{representatives}
      """
}

process pMapBwa {
    label 'large'
    container "${params.samtools_bwa_image}"
    when params.steps.containsKey("readMapping")
    publishDir params.output, saveAs: { filename -> getOutput("${bin_shuffle_id}",params.runid ,"bwa", filename) }
    input:
      tuple path(sample), val(bin_shuffle_id), path(representatives_fasta), path(x, stageAs: "*") 
    output:
      tuple val("${ID}"), val(bin_shuffle_id), path("*bam"), path("*bam.bai"), emit: alignment
      tuple file(".command.sh"), file(".command.out"), file(".command.err"), file(".command.log")
    shell:
    template('bwa.sh')
}


process pCovermCount {
    when params.steps.containsKey("readMapping")
    label 'small'
    publishDir params.output, saveAs: { filename -> getOutput("${sample}", params.runid, "coverm", filename) }
    input:
      tuple val(sample), file(mapping), file(index), file(list_of_representatives)
    output:
      tuple val("${sample}"), path("${sample}_out/mean.tsv"), emit: mean
      tuple val("${sample}"), path("${sample}_out/trimmed_mean.tsv"), emit: trimmedMean
      tuple val("${sample}"), path("${sample}_out/count.tsv"), emit: count
      tuple val("${sample}"), path("${sample}_out/rpkm.tsv"), emit: rpkm
      tuple val("${sample}"), path("${sample}_out/tpm.tsv"), emit: tpm
      tuple val("${sample}"), path("${sample}_out/mean_mincov10.tsv"), emit: meanMincov10
      tuple file(".command.sh"), file(".command.out"), file(".command.err"), file(".command.log")
    shell:
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
       | set {samples}

     _wReadMappingBwa(samples, genomesList)
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
     samples
     genomes
   main:
     _wReadMappingBwa(samples, genomes)
   emit:
     trimmedMean = _wReadMappingBwa.out.trimmedMean
}


workflow _wReadMappingBwa {
   take:
     samples
     genomes
   main:
     // Create temporary directory for merged fasta files
     genomesTempDir = params.tempdir + "/genomesToBeMerged"
     file(genomesTempDir).mkdirs()

     // Concatenate genomes 
     genomes | collectFile(tempDir: genomesTempDir){ item -> [ "mergedGenomes.fasta", item.text ] } \
      | set { genomesMerged }

     // Create BWA index of all genomes
     SAMPLE_IDX=1
     BWA_INDEX_IDX=0
     GENOMES_IDX=2
     genomesMerged | pBwaIndex | map{ bwaIndex -> [bwaIndex]} | combine(samples) \
      | combine(genomesMerged) \
      | map{ it -> [it[SAMPLE_IDX].READS, it[SAMPLE_IDX].SAMPLE, it[GENOMES_IDX], it[BWA_INDEX_IDX]] } \
      | set {index}

     pMapBwa(index)
     // Map all samples against all genomes using bwa 
     pMapBwa.out.alignment | combine(genomes | map {it -> file(it)} \
      | toList() | map { it -> [it]}) | pCovermCount
   emit:
     trimmedMean = pCovermCount.out.trimmedMean
}
