nextflow.enable.dsl=2

process runMetaSpades {

    label 'large'

    tag "$sample"

    container 'quay.io/biocontainers/spades:${metaspades_tag}'

    publishDir "${params.output}/${sample}/assembly/metaspades/${metaspades_tag}" 

    when params.metaspades

    input:
    tuple val(sample), path(fastqs, stageAs: 'fastq.fq.gz')

    output:
    tuple val("${sample}"), env(TYPE), path("${sample}/contigs.fasta"), emit: contigs

    shell:
    '''
    metaspades.py --threads 28  --12 fastq.fq.gz -o !{sample}
    TYPE="metaspades"
    '''
}


process runMegahit {

    label 'large'

    tag "$sample"

    container "vout/megahit:${megahit_tag}"

    publishDir "${params.output}/${sample}/assembly/megahit/${megahit_tag}" 

    when params.containsKey("megahit") && params.assembly == "megahit"

    cpus 28
    input:
    tuple val(sample), path(fastqs, stageAs: 'fastq.fq.gz')

    output:
    tuple val("${sample}"), env(TYPE), path("${sample}/final.contigs.fa"), emit: contigs

    shell:
    '''
    megahit -t !{task.cpus} --12 !{fastqs}
    TYPE="megahit" 
    mkdir !{sample}
    mv megahit_out/final.contigs.fa !{sample}/
    '''
}


process runTrimmomatic {

    label 'large'

    tag "$sample"

    publishDir "${params.output}/${sample}/reads/" 

    container "quay.io/biocontainers/trimmomatic:${params.trimmomatic_tag}"

    input:
    tuple val(sample), path(genomeReads1), path(genomeReads2)

    output:
    tuple val("${sample}"), path("*_R1.p.fastq.gz"), path("*_R2.p.fastq.gz")

    script:

    fq_1_paired = sample + '_R1.p.fastq.gz'
    fq_1_unpaired = sample + '_R1.s.fastq.gz'
    fq_2_paired = sample + '_R2.p.fastq.gz'
    fq_2_unpaired = sample + '_R2.s.fastq.gz'

    """
    trimmomatic \
    PE -phred33 \
    -threads ${task.cpus} \
    ${genomeReads1} \
    ${genomeReads2} \
    $fq_1_paired \
    $fq_1_unpaired \
    $fq_2_paired \
    $fq_2_unpaired \
    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36
    """
}



process runFastp {

    label 'large'

    tag "$sample"

    publishDir "${params.output}/${sample}/reads/fastp/${params.fastp_tag}"

    container "quay.io/biocontainers/fastp:${params.fastp_tag}"

    input:
    tuple val(sample), path(genomeReads1), path(genomeReads2)

    output:
    tuple val("${sample}"), path("*_R1.p.fastq.gz"), path("*_R2.p.fastq.gz"), emit: fastq
    path("*_fastp_report.html"), emit: html_report
    path("*_fastp_report.json"), emit: json_report
    path("fastp_summary"), emit: fastp_summary
    

    script:
    fq_1_paired = sample + '_R1.p.fastq.gz'
    fq_1_unpaired = sample + '_R1.s.fastq.gz'
    fq_2_paired = sample + '_R2.p.fastq.gz'
    fq_2_unpaired = sample + '_R2.s.fastq.gz'
    template 'fastp.sh'
}


process runBBMapInterleave {

    label 'large'

    tag "$sample"

    publishDir "${params.output}/${sample}/reads/bbmap/" 

    container "quay.io/biocontainers/bbmap:${params.bbmap_tag}"

    input:
    tuple val(sample), path(read1, stageAs: "read1.fq.gz"), path(read2, stageAs: "read2.fq.gz")

    output:
    tuple val("${sample}"), path("reads.fq.gz")

    shell:
    """
    reformat.sh in1=read1.fq.gz in2=read2.fq.gz out=reads.fq.gz
    """
}

process runMegahitInterleaved {

    label 'large'

    tag "$sample"

    publishDir "${params.output}/${sample}/assembly/megahit/${params.megahit_tag}" 

    when params.containsKey("megahit") && params.assembly == "megahit"

    errorStrategy 'ignore'

    container "vout/megahit:release-v1.2.9:${params.megahit_tag}"

    input:
    tuple val(sample), path(interleaved_reads, stageAs: "interleaved.fq.gz")

    output:
    tuple val("${sample}"), env(TYPE), path("final.contigs.fa"), emit: contigs
    tuple val("${sample}"), path("interleaved.fastp.fq.gz"), emit: reads_processed
    tuple val("${sample}"), path("fastp_summary.tsv"), emit: fastp_summary


    shell:
    template 'megahit_interleaved.sh'

}


process runMegahitSplit {

    label 'large'

    tag "$sample"

    publishDir "${params.output}/${sample}/assembly/megahit/${params.megahit_tag}" 

    errorStrategy 'ignore'

    when params.containsKey("megahit") && params.assembly == "megahit"

    container "vout/megahit:${params.megahit_tag}"

    input:
    tuple val(sample), path(read1, stageAs: "read1.fq.gz"), path(read2, stageAs: "read2.fq.gz")

    output:
    tuple val("${sample}"), env(TYPE), path("final.contigs.fa"), emit: contigs
    tuple val("${sample}"), path("interleaved.fastp.fq.gz"), emit: reads_processed
    tuple val("${sample}"), path("fastp_summary_before.tsv"), emit: fastp_summary_before
    tuple val("${sample}"), path("fastp_summary_after.tsv"), emit: fastp_summary_after
    tuple val("${sample}"), path("fastp.json"), emit: fastp_summary

    shell:
    template 'megahit_split.sh'

}

process runBBMapDeinterleave {

    label 'large'

    tag "$sample"

    publishDir "${params.output}/${sample}/reads/bbmap/${params.bbmap_tag}" 

    container "quay.io/biocontainers/bbmap:${params.bbmap_tag}"

    input:
    tuple val(sample), path(interleaved_reads, stageAs: "interleaved.fq.gz")

    output:
    tuple val("${sample}"), path("read1.fq.gz"), path("read2.fq.gz")

    shell:
    """
    reformat.sh in=interleaved.fq.gz out1=read1.fq.gz out2=read2.fq.gz
    """
}


workflow run_assembly {
   take: 
     input_split_raw_reads
   main:
     input_split_raw_reads | runFastp 
     runFastp.out.fastq | runBBMapInterleave | set{input_reads}

     runFastp.out.fastp_summary_before | collectFile(newLine: false, keepHeader: true, storeDir: params.output + "/summary/" ){ item ->
          [ "fastp_summary_before.tsv", item[1].text  ]
     }

     runFastp.out.fastp_summary_after | collectFile(newLine: false, keepHeader: true, storeDir: params.output + "/summary/" ){ item ->
          [ "fastp_summary_after.tsv", item[1].text  ]
     }
 
      
     input_reads | runMegahit | set { megahit }
     input_reads | runMetaSpades | set { metaspades }

     metaspades.contigs | mix(megahit.contigs) |  set{contigs}
   emit:
     contigs = contigs
     processed_reads = input_reads
     fastp_summary = runFastp.out.fastp_summary
}

