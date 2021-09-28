nextflow.enable.dsl=2

params.samples = "/raid/CAMI/CAMI_MOUSEGUT/readmapping/workflow_all_reads/samples.tsv"
params.mapping = "/vol/spool/CAMI/CAMI_MOUSEGUT/mapping/workflow/representatives_sample1.tsv_mapping"
params.representatives = "/vol/spool/CAMI/CAMI_MOUSEGUT/mapping/workflow/representatives_sample1_test.fa"  
params.number_samples = 20
params.all_representatives = "/raid/CAMI/CAMI_MOUSEGUT/readmapping/workflow_all_reads/representatives.tsv" 
params.shuffle = 400

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
      tuple val(id), path(representatives)
    output:
      tuple val("${id}") ,path('*.{amb,ann,bwt,pac,sa,fa}')
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
      tuple val(mode), path(sample), val(bin_shuffle_id), val(ID), path(representatives_fasta), path(x, stageAs: "*") 
    output:
      tuple val("${ID}"), val(bin_shuffle_id), path("*bam"), emit: alignmentIndex
      tuple val("${ID}"), val(bin_shuffle_id), path("*bam"), emit: alignment
      tuple file(".command.sh"), file(".command.out"), file(".command.err"), file(".command.log")
    shell:
    MODE = mode == "paired" ? " -p " : "" 
    template('bwa.sh')
}

process pMapBwaCami {
    label 'large'
    container "${params.samtools_bwa_image}"
    publishDir params.output, saveAs: { filename -> getOutput("${sample}", params.runid ,"bwa" , filename) }
    input:
      tuple val(mode), path(sample), val(bin_shuffle_id), val(ID), path(representatives_fasta), path(x, stageAs: "*") 
    output:
      tuple val("${ID}"), val(bin_shuffle_id), path("*.sam.gz")
      tuple file(".command.sh"), file(".command.out"), file(".command.err"), file(".command.log")
    shell:
    template('bwa.sh')
}

process pBwaCount {
    label 'tiny'
    publishDir params.output, saveAs: { filename -> getOutput("${sample}", params.runid, "coverm", filename) }
    input:
      tuple val(ID), file(sample), file(mapping)
    output:
      path "${sample}.count"
      tuple file(".command.sh"), file(".command.out"), file(".command.err"), file(".command.log")
    shell:
    template('count.sh')
}

process pCovermCount {
    when params.steps.containsKey("readMapping")
    label 'small'
    publishDir params.output, saveAs: { filename -> getOutput("${sample}", params.runid, "coverm", filename) }
    input:
      tuple val(sample), file(mapping), file(index), file(list_of_representatives)
    output:
      path("${sample}_out", type: "dir")
      tuple file(".command.sh"), file(".command.out"), file(".command.err"), file(".command.log")
    shell:
    template('coverm.sh')
}


process pMergeAlignment {
    when params.steps.containsKey("readMapping")
    label 'tiny'
    container 'quay.io/biocontainers/samtools:1.12--h9aed4be_1'
    publishDir params.output, saveAs: { filename -> getOutput("${sample}", params.runid, "coverm", filename) }
    input:
      tuple val(ID), val(sample), file("alignment?.bam")
    output:
      tuple val("${sample}"), path("${sample}.bam"), path("*bam.bai"), emit: alignmentIndex
      tuple file(".command.sh"), file(".command.out"), file(".command.err"), file(".command.log")
    shell:
    """
    samtools merge !{sample}.bam alignment*.bam
    samtools index !{sample}.bam
    """
}

process pShuff {
    publishDir params.output + "/sample"
    label 'tiny'
    input:
      path file_of_bins
      each samples_number
    output:
      path "${file_of_bins}_${samples_number}_shuffled.tsv", emit: tsv
    shell:
    '''
    shuf -n !{params.shuffle} !{file_of_bins} > !{file_of_bins}_!{samples_number}_shuffled.tsv
    '''
}

process pCreateMapping {
    publishDir params.output + "/mapping"
    label 'tiny'
    input:
      path representatives_sample
    output:
      tuple val("${representatives_sample}"), path("${representatives_sample}_mapping"), emit: mapping_out
      tuple val("${representatives_sample}"), path("${representatives_sample}.fa"), emit: fa
    shell:
    '''
    mkdir -p representatives; cat  !{representatives_sample} | xargs -P 20 -I {}  sh -c 'bin="$1";bin_name=$(basename $bin); grep ">" $bin | tr -d ">" | sed "s/$/\t$bin_name/g" > representatives/$bin_name ' sh {} ; cat representatives/* > !{representatives_sample}_mapping
    cat !{representatives_sample} | xargs -I {} cat {} > !{representatives_sample}.fa
    '''
}

workflow pCreateSamples {
   main:
     shuff(Channel.value(params.all_representatives), 1..params.number_samples)
     create_mapping(shuff.out.tsv)
   emit:
     bin_sample_mapping_fa = create_mapping.out.fa
     bin_sample_mapping_tsv = create_mapping.out.mapping_out
}


workflow wBwaMultiple {
    take:
      representatives_file
      samples
      mapping
    main:
      samples \
       | splitCsv(sep: '\t', header: true) | set {samples_split}
      representatives_file | pBwaIndex | combine(samples_split) \
       | combine(representatives_file, by: 0) \
       | map{ it -> [it[2].READS, it[2].SAMPLE, it[0], it[3], it[1]] } \
       | set { index }
      pMapBwaCami(index) | set { sam_gz  }
      sam_gz | combine(mapping, by: 0) | map { it -> it[(1..3)] } | pBwaCount

}

workflow wReadMappingBwa {
   take:
     id
     representatives
     samplesPaired
     samplesSingle
     representativesList
   main:
     samplesPaired | splitCsv(sep: '\t', header: true) | set {samplesPairedSplit}

     samplesSingle | splitCsv(sep: '\t', header: true) | set {samplesSingleSplit}

     id | combine(representatives) | pBwaIndex | set {index} 

     index | combine(samplesPairedSplit) \
      | combine(representatives) \
      | map{ it -> ["paired", it[2].READS, it[2].SAMPLE, it[0], it[3], it[1]] } \
      | set {paired}

    index | combine(samplesSingleSplit) \
      | combine(representatives) \
      | map{ it -> ["single", it[2].READS, it[2].SAMPLE, it[0], it[3], it[1]] } \
      | set {single}

     // Paired and single reads should be mapped back
     pMapBwa(paired | mix(single))
     
     // The resulting alignments (bam files) should merged
     pMapBwa.out.alignment | groupTuple(by: 1) | pMergeAlignment
         
     // The alignment is then analysed per dataset
     pMergeAlignment.out.alignmentIndex | combine(representativesList \
      | map {it -> file(it)} | toList() | map { it -> [it]}) | pCovermCount
}

workflow {
   wCreateSamples() 
   wBwaMultiple(wCreateSamples.out.bin_sample_mapping_fa, channel.fromPath(params.samples), create_samples.out.bin_sample_mapping_tsv)
}
