nextflow.enable.dsl=2

params.samples = "/raid/CAMI/CAMI_MOUSEGUT/readmapping/workflow_all_reads/samples.tsv"
params.out = "out"
params.mapping = "/vol/spool/CAMI/CAMI_MOUSEGUT/mapping/workflow/representatives_sample1.tsv_mapping"
params.representatives = "/vol/spool/CAMI/CAMI_MOUSEGUT/mapping/workflow/representatives_sample1_test.fa"  
params.number_samples = 20
params.all_representatives = "/raid/CAMI/CAMI_MOUSEGUT/readmapping/workflow_all_reads/representatives.tsv" 
params.shuffle = 400

process bwa_index {
    conda ' bioconda::bwa=0.7.17'
    label 'large'
    input:
      tuple val(id), path(representatives)
    output:
      tuple val("${id}") ,path('*.{amb,ann,bwt,pac,sa,fa}')
    shell:
      """
      bwa index !{representatives}
      """
}

process map_bwa {
    label 'large'
    conda 'bioconda::bwa=0.7.17'
    publishDir params.out + "/bwa"
    input:
      tuple path(sample), val(bin_shuffle_id), val(ID), path(representatives_fasta), path(x, stageAs: "*") 
    output:
      tuple val("${ID}"), val(bin_shuffle_id), path("*bam.sorted"), path("*bam.sorted.bai")
    shell:
    template('bwa.sh')
}

process map_bwa_cami {
    label 'large'
    conda 'bioconda::bwa=0.7.17'
    publishDir params.out + "/bwa"
    input:
      tuple path(sample), val(bin_shuffle_id), val(ID), path(representatives_fasta), path(x, stageAs: "*") 
    output:
      tuple val("${ID}"), val(bin_shuffle_id), path("*.sam.gz")
    shell:
    template('bwa.sh')
}

process bwa_count {
    label 'tiny'
    publishDir params.out + "/count"
    input:
      tuple val(ID), path(sample), path(mapping)
    output:
      path "${sample}.count"
    shell:
    template('count.sh')
}

process coverm_count {
//   conda 'bioconda::coverm'
    label 'small'
    publishDir params.out + "/count"
    input:
      tuple val(ID), val(sample), path(mapping), path(index), path(list_of_representatives)
    output:
      path("${ID}_${sample}_out", type: "dir")
    shell:
    template('coverm.sh')
}

process shuff {
    publishDir params.out + "/sample"
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

process create_mapping {
    publishDir params.out + "/mapping"
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

workflow create_samples {
   main:
     shuff(Channel.value(params.all_representatives), 1..params.number_samples)
     create_mapping(shuff.out.tsv)
   emit:
     bin_sample_mapping_fa = create_mapping.out.fa
     bin_sample_mapping_tsv = create_mapping.out.mapping_out
}


workflow bwa_multiple {
    take:
      representatives_file
      samples
      mapping
    main:
      samples \
       | splitCsv(sep: '\t', header: true) | set {samples_split}
      representatives_file | bwa_index | combine(samples_split) \
       | combine(representatives_file, by: 0) \
       | map{ it -> [it[2].READS, it[2].SAMPLE, it[0], it[3], it[1]] } \
       | view() \
       | set { index }
      map_bwa_cami(index) | set { sam_gz  }
      sam_gz | combine(mapping, by: 0) | map { it -> it[(1..3)] } | bwa_count

}

workflow bwa {
   take: 
     id
     representatives
     samples
     list_of_representatives
   main:
     samples | splitCsv(sep: '\t', header: true) | set {samples_split}
     id | combine(representatives) | bwa_index | combine(samples_split) \
      | combine(representatives) \
      | map{ it -> [it[2].READS, it[2].SAMPLE, it[0], it[3], it[1]] } \
      | set {index}
     map_bwa(index) | combine(list_of_representatives) | view() | coverm_count

}

workflow {
   create_samples() 
   bwa_multiple(create_samples.out.bin_sample_mapping_fa, channel.fromPath(params.samples), create_samples.out.bin_sample_mapping_tsv)
}
