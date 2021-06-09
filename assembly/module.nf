nextflow.enable.dsl=2

process pMetaSpades {

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


process pMegahit {

    label 'large'

    tag "$sample"

    container "vout/megahit:${megahit_tag}"

    publishDir "${params.output}/${sample}/assembly/megahit/${megahit_tag}" 

    when params.steps.assembly.containsKey("megahit")

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


process pTrimmomatic {

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



process pFastp {

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


process pBBMapInterleave {

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

process pMegahitInterleaved {

    label 'large'

    tag "$sample"

    publishDir "${params.output}/${sample}/assembly/megahit/${params.megahit_tag}" 

    when params.steps.mags_generation.assembly == "megahit"

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


process pMegahitSplit {

    label 'large'

    tag "$sample"

    publishDir "${params.output}/${sample}/assembly/megahit/${params.megahit_tag}" 

    errorStrategy 'ignore'

    when params.steps.assembly.containsKey("megahit")

    container "vout/megahit:${params.megahit_tag}"

    input:
    tuple val(sample), path(read1, stageAs: "read1.fq.gz"), path(read2, stageAs: "read2.fq.gz")

    output:
    tuple val("${sample}"), env(TYPE), path("${sample}_final.contigs.fa.gz"), emit: contigs
    tuple val("${sample}"), path("interleaved.fastp.fq.gz"), emit: reads_processed
    tuple val("${sample}"), path("fastp_summary_before.tsv"), emit: fastp_summary_before
    tuple val("${sample}"), path("fastp_summary_after.tsv"), emit: fastp_summary_after
    tuple val("${sample}"), path("fastp.json"), emit: fastp_summary
    tuple val("${sample}"), path("*_report.html"), emit: fastp_summary_html
    tuple val("${sample}"), path("${sample}_contigs_stats.tsv"), emit: contigs_stats

    shell:
    template 'megahit_split.sh'
}

process pBBMapDeinterleave {

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


workflow wAssembly {
   take: 
     input_split_raw_reads
   main:
     input_split_raw_reads | pFastp 
     pFastp.out.fastq | pBBMapInterleave | set{input_reads}

     pFastp.out.fastp_summary_before | collectFile(newLine: false, keepHeader: true, storeDir: params.output + "/summary/" ){ item ->
          [ "fastp_summary_before.tsv", item[1].text  ]
     }

     pFastp.out.fastp_summary_after | collectFile(newLine: false, keepHeader: true, storeDir: params.output + "/summary/" ){ item ->
          [ "fastp_summary_after.tsv", item[1].text  ]
     }
 
      
     input_reads | pMegahit | set { megahit }
     input_reads | pMetaSpades | set { metaspades }

     metaspades.contigs | mix(megahit.contigs) |  set{contigs}
   emit:
     contigs = contigs
     processed_reads = input_reads
     fastp_summary = runFastp.out.fastp_summary
}


// Fastq interleaved, QC and Assembler are executed in one process
workflow _wFastqSplitSkip {
       take:
         reads_table
       main:
            reads_table | splitCsv(sep: '\t', header: true) \
             | map { it -> [ it.SAMPLE, it.READS1, it.READS2 ]} \
             | pMegahitSplit 
             
            pMegahitSplit.out.reads_processed | set { processed_reads}

            pMegahitSplit.out.fastp_summary | collectFile(newLine: false, keepHeader: true, storeDir: params.output + "/summary/" ){ item ->
              [ "fastp_summary.tsv", item[1].text ]
            }

            pMegahitSplit.out.fastp_summary_after | collectFile(newLine: false, keepHeader: true, storeDir: params.output + "/summary/" ){ item ->
              [ "fastp_summary_after.tsv", item[1].text ]
            }

            pMegahitSplit.out.fastp_summary_before | collectFile(newLine: false, keepHeader: true, storeDir: params.output + "/summary/" ){ item ->
              [ "fastp_summary_before.tsv", item[1].text ]
            }

            pMegahitSplit.out.contigs_stats | collectFile(newLine: false, keepHeader: true, storeDir: params.output + "/summary/" ){ item ->
              [ "contigs_stats.tsv", item[1].text ]
            }

      emit:
        contigs = pMegahitSplit.out.contigs
        processed_reads = processed_reads
}


// Fastq split, QC and Assembler are executed in sepearate processes
workflow _wFastqSplitSeperate {
       take:
         reads_table
       main:
         reads_table | splitCsv(sep: '\t', header: true) \
          | map { it -> [ it.SAMPLE, it.READS1, it.READS2 ]} \
          | wAssembly

         run_assembly.out.fastp_summary | collectFile(newLine: false, keepHeader: true, storeDir: params.output + "/summary/"  ){ item ->
            [ "fastp_summary.tsv", item[1].text  ]
         }

} 

// Fastq interleaved, QC and Assembler are executed in once process
workflow _wFastqInterleavedSkip {
       take:
         reads_table
       main:
         reads_table | splitCsv(sep: '\t', header: true) \
          | map { it -> [ it.SAMPLE, it.READS ]} \
          | pMegahitInterleaved     

         pMegahitInterleaved.out.reads_processed | set { processed_reads}
         pMegahitInterleaved.out.fastp_summary | collectFile(newLine: false, keepHeader: true, storeDir: params.output + "/summary/" ){ item ->
           [ "fastp_summary.tsv", item[1].text  ]
         }
}


// Fastq interleaved, QC and Assembler are executed in seperate processes
workflow _wFastqInterleavedSeperate {
       take:
         reads_table
       main:
         reads_table | splitCsv(sep: '\t', header: true) \
          | map { it -> [ it.SAMPLE, it.READS ]} \
          | pBBMapDeinterleave | pAssembly
         pAssembly.out.fastp_summary | collectFile(newLine: false, keepHeader: true, storeDir: params.output + "/summary/"  ){ item ->
           [ "fastp_summary.tsv", item[1].text  ]
         }
}


/*
 * Takes a tab separated file of files containing reads as input and produces assembly, 
 * binning results and contamination, completeness and taxonomy reports for mags.
 * Input file with columns seperated by tabs:
 * Dataset_ID Left_Read Right_Read
 *
 * Left and right read could be https, s3 links or file path.
 * 
 */
workflow wAssemblyFile {
     take:
       reads
     main:
       if(params.steps.assembly.interleaved){
          if(params.steps.assembly.skip){
            _wFastqInterleavedSkip(reads)
          } else {
            _wFastqInterleavedSeperate(reads)
          }
       }

       if(!params.steps.assembly.interleaved){
          if(params.steps.assembly.skip){
           results = _wFastqSplitSkip(reads)
          } else {
           _wFastqSplitSeperate(reads)
          }
      }
    emit:
      contigs = results.contigs
      processed_reads = results.processed_reads
}
