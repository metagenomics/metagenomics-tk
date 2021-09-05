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
    container "quay.io/biocontainers/bwa:${params.bwa_tag}"
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
    container "pbelmann/bwa-samtools:${params.samtools_bwa_tag}"
    when params.steps.containsKey("readMapping")
    publishDir params.output, saveAs: { filename -> getOutput("${bin_shuffle_id}",params.runid ,"bwa", filename) }
    input:
      tuple path(sample), val(bin_shuffle_id), val(ID), path(representatives_fasta), path(x, stageAs: "*") 
    output:
      tuple val("${ID}"), val(bin_shuffle_id), path("*bam"), path("*bam.bai")
    shell:
    template('bwa.sh')
}

process pMapBwaCami {
    label 'large'
    container "pbelmann/bwa-samtools:${params.samtools_bwa_tag}"
    publishDir params.output, saveAs: { filename -> getOutput("${sample}", params.runid ,"bwa" , filename) }
    input:
      tuple path(sample), val(bin_shuffle_id), val(ID), path(representatives_fasta), path(x, stageAs: "*") 
    output:
      tuple val("${ID}"), val(bin_shuffle_id), path("*.sam.gz")
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
    shell:
    template('count.sh')
}

process pCovermCount {
    when params.steps.containsKey("readMapping")
    label 'small'
    publishDir params.output, saveAs: { filename -> getOutput("${sample}", params.runid, "coverm", filename) }
    input:
      tuple val(ID), val(sample), file(mapping), file(index), file(list_of_representatives)
    output:
      path("${ID}_${sample}_out", type: "dir")
    shell:
    template('coverm.sh')
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
     samples
     representativesList
   main:
     samples | splitCsv(sep: '\t', header: true) | set {samples_split}
     id | combine(representatives) | pBwaIndex | combine(samples_split) \
      | combine(representatives) \
      | map{ it -> [it[2].READS, it[2].SAMPLE, it[0], it[3], it[1]] } \
      | set {index}
     pMapBwa(index) | combine(representativesList | map {it -> file(it)} | toList() | map { it -> [it]}) | pCovermCount
}

workflow {
   wCreateSamples() 
   wBwaMultiple(wCreateSamples.out.bin_sample_mapping_fa, channel.fromPath(params.samples), create_samples.out.bin_sample_mapping_tsv)
}
