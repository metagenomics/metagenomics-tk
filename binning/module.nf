nextflow.enable.dsl=2
params.maxbin = false


process pGetMappingQuality {

    container 'quay.io/biocontainers/samtools:1.12--h9aed4be_1'

    tag "$sample"

    publishDir "${params.output}/${sample}/mapping/bowtie/"

    errorStrategy 'retry'

    label 'tiny'

    input:
    tuple val(sample), val(TYPE), path(bam)

    output:
    tuple val("${sample}"), val("${TYPE}"), file("${sample}_${TYPE}_flagstat.tsv"), emit: flagstat_raw
    tuple val("${sample}"), val("${TYPE}"), file("${sample}_${TYPE}_flagstat_passed.tsv"), emit: flagstat_passed
    tuple val("${sample}"), val("${TYPE}"), file("${sample}_${TYPE}_flagstat_failed.tsv"), emit: flagstat_failed

    shell:
    template 'mapping_quality.sh'
}



process pBowtie {

    container "pbelmann/bowtie2:${params.bowtie_tag}"

    label 'large'

    tag "$sample"

    publishDir "${params.output}/${sample}/mapping/bowtie/${params.bowtie_tag}" 

    errorStrategy 'retry'

    cpus 28

    input:
    tuple val(sample), val(TYPE), path(contigs), path(fastqs, stageAs: 'fastq.fq.gz')

    output:
    tuple val("${sample}"), val("${TYPE}"), file("${TYPE}_${sample}.bam"), emit: bam

    shell:
    '''
    INDEX=!{sample}.index
    bowtie2-build --threads 28 --quiet !{contigs} $INDEX 
    bowtie2 -p 28 --quiet --very-sensitive -x $INDEX --interleaved fastq.fq.gz | samtools view -F 3584 --threads 28 -bS - | samtools sort -l 9 --threads 28 - > !{TYPE}_!{sample}.bam
    '''
}


process pMetabat {

    container "metabat/metabat:${params.metabat_tag}"

    errorStrategy 'ignore'

    tag "$sample"

    label 'large'

    publishDir "${params.output}/${sample}/binning/metabat/${params.metabat_tag}" 

    when params.steps.binning.containsKey("metabat")

    input:
    tuple val(sample), val(TYPE), path(contigs), path(bam)

    output:
    tuple val("${sample}"), env(NEW_TYPE), file("${TYPE}/*bin*"), optional: true, emit: bins
    tuple val("${sample}"), env(NEW_TYPE), file("*.depth.txt"), optional: true, emit: metabat_depth
    tuple val("${sample}"), env(NEW_TYPE), file("${sample}_bins_depth.tsv"), optional: true, emit: bins_depth
    tuple val("${sample}"), env(NEW_TYPE), file("${sample}_bins_stats.tsv"), optional: true, emit: bins_stats

    shell:
    template 'metabat.sh'
}


process pMaxBin {

    container "quay.io/biocontainers/maxbin2:${params.maxbin_tag}"

  //  errorStrategy 'ignore'

    label 'large'

    tag "$sample"

    when params.maxbin

    publishDir "${params.output}/${sample}/binning/maxbin/${params.maxbin_tag}" 

    input:
    tuple val(sample), val(TYPE), path(contigs), path(reads)

    output:
    tuple val("${sample}"), env(NEW_TYPE), file("${TYPE}_${sample}/out.*.fasta"), optional: true, emit: bins

    shell:
    '''
    NEW_TYPE="!{TYPE}_maxbin"
    mkdir !{TYPE}_!{sample}
    run_MaxBin.pl -preserve_intermediate -contig !{contigs} -reads !{reads} -thread 28 -out !{TYPE}_!{sample}/out
    '''
}


def bufferBins(binning){
  def chunkList = [];
  binning[2].each {
     chunkList.add([binning[0],binning[1], it]);
  }
  return chunkList;
}


def aslist(element){
  if(element instanceof Collection){
    return element;
  } else {
    return [element];
  }
}


workflow wBinning {
   take: 
     contigs
     input_reads
   main:
     contigs | join(input_reads | mix(input_reads), by: 0) | pBowtie | set { bam }
     bam | pGetMappingQuality 
     
     pGetMappingQuality.out.flagstat_passed | collectFile(newLine: false, keepHeader: true, storeDir: params.output + "/summary/"){ item ->
       [ "flagstat_passed.tsv", item[2].text  ]
     }

     pGetMappingQuality.out.flagstat_failed | collectFile(newLine: false, keepHeader: true, storeDir: params.output + "/summary/"){ item ->
       [ "flagstat_failed.tsv", item[2].text  ]
     }

     contigs | join(bam, by: [0, 1]) | pMetabat | set { metabat }
     contigs | join(input_reads | mix(input_reads), by: 0) | pMaxBin | set { maxbin }

     metabat.bins | map({ it -> it[2] = aslist(it[2]); it  }) | set{ bins_list }
     metabat.bins_stats | map { it -> file(it[2]) } | splitCsv(sep: '\t', header: true) | set { bins_stats }
     metabat.bins_depth | collectFile(newLine: false, keepHeader: true, storeDir: params.output + "/summary/"){ item ->
       [ "metabat_bins_depth.tsv", item[2].text  ]
     }

     metabat.bins_stats | collectFile(newLine: false, keepHeader: true, storeDir: params.output + "/summary/"){ item ->
       [ "metabat_bins_depth.tsv", item[2].text  ]
     }

    // metabat.bins | mix(maxbin.bins) | set {bins}
   emit:
     bins_stats = bins_stats
     bins = bins_list
     mapping = bam
}


