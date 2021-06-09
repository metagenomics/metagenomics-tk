nextflow.enable.dsl=2
include { wUnmappedReadsFile } from '../../sampleAnalysis/module'


process pFrHit {

    errorStrategy 'ignore'

    label 'large'

    tag "$sample"

    stageInMode 'copy'

    publishDir "${params.output}/${sample}/fragmentRecruitment"

    container "pbelmann/bwa-samtools:${params.samtools_bwa_tag}"

    when params.steps.containsKey("fragmentRecruitment")

    input:
    tuple val(sample), file(reads), file(genomesCombined), file(genomesList) 

    output:
    tuple val("${sample}"), file("${sample}.bam"), optional: true, emit: alignment
    tuple val("${sample}"), file("coverage"), optional: true, emit: coverageStats

    shell:
    '''
    mkdir coverage
    seqkit fq2fa !{reads} -o reads.fasta.gz
    AVG_LEN=$(seqkit stats -T reads.fasta.gz | tail -n 1 | cut -f 7)
    HALF_AVG_LEN=$(bc <<< "${AVG_LEN} / 2")
    gunzip reads.fasta.gz
    echo "Min. Coverage Parameter: ${HALF_AVG_LEN}"
    fr-hit -c 95 -f 1 -m ${HALF_AVG_LEN} -T !{task.cpus} -a reads.fasta -d !{genomesCombined} -o coverage/out.psl
    if [ -s "coverage/out.psl" ] 
    then
      psl2sam.pl coverage/out.psl | samtools view -bT !{genomesCombined} - | samtools calmd -E - genomes | samtools sort > !{sample}.bam 2> /dev/null
      coverm genome -t !{task.cpus} --min-covered-fraction 0 -b !{sample}.bam --genome-fasta-list !{genomesList} --methods count --output-file coverage/readCount.tsv
      coverm genome -t !{task.cpus} --min-covered-fraction 0 -b !{sample}.bam --genome-fasta-list !{genomesList} --methods covered_fraction --output-file coverage/coveredFraction.tsv
      coverm genome -t !{task.cpus} --min-covered-fraction 0 -b !{sample}.bam --genome-fasta-list !{genomesList} --methods covered_bases --output-file coverage/coveredBases.tsv
      coverm genome -t !{task.cpus} --min-covered-fraction 0 -b !{sample}.bam --genome-fasta-list !{genomesList} --methods length --output-file coverage/genomeLength.tsv
      coverm genome -t !{task.cpus} --min-covered-fraction 0 -b !{sample}.bam --genome-fasta-list !{genomesList} --methods trimmed_mean --output-file coverage/trimmedMean.tsv
    else
      echo "No reads could be recruited!"
    fi
    '''
}



process pUnzip {

  label 'tiny'

  input:
  file x

  output:
  file("${x.baseName}")

  script:
  """
  < $x zcat > ${x.baseName}
  """
}



process pCombinedAlignmentAnalysis {

    errorStrategy 'ignore'

    label 'medium'

    when params.steps.containsKey("fragmentRecruitment")

    publishDir "${params.output}/fragmentRecruitment"

    container "pbelmann/bwa-samtools:${params.samtools_bwa_tag}"

    input:
    file(alignments)
    tuple file(genomesCombined), file(genomesList)

    output:
    path("coverage"), emit: coverageStats

    shell:
    '''
    mkdir coverage
    samtools merge -@ !{task.cpus} combined_alignments.bam !{alignments}
    coverm genome -t !{task.cpus} --min-covered-fraction 0 -b combined_alignments.bam --genome-fasta-list !{genomesList} --methods count --output-file coverage/readCount.tsv
    coverm genome -t !{task.cpus} --min-covered-fraction 0 -b combined_alignments.bam --genome-fasta-list !{genomesList} --methods covered_fraction --output-file coverage/coveredFraction.tsv
    coverm genome -t !{task.cpus} --min-covered-fraction 0 -b combined_alignments.bam --genome-fasta-list !{genomesList} --methods covered_bases --output-file coverage/coveredBases.tsv
    coverm genome -t !{task.cpus} --min-covered-fraction 0 -b combined_alignments.bam --genome-fasta-list !{genomesList} --methods length --output-file coverage/genomeLength.tsv
    coverm genome -t !{task.cpus} --min-covered-fraction 0 -b combined_alignments.bam --genome-fasta-list !{genomesList} --methods trimmed_mean --output-file coverage/trimmedMean.tsv
    '''

}


workflow wFragmentRecruitmentFile {
   take: 
     sampleReadsFile
     genomes
   main:
     fragmentRecruitmentGenomes = params.tempdir + "/fragmentRecruitmentGenomes"
     file(fragmentRecruitmentGenomes).mkdirs()

     sampleReadsFile | splitCsv(sep: '\t', header: true) \
       | map { sample -> [sample.SAMPLE, file(sample.READS)] } | set {sampleReadsList}

     _wFragmentRecruitment(sampleReadsList, genomes)
}


workflow wFragmentRecruitmentList {
   take:
     sampleReads
     genomes
   main:
     SAMPLE_NAME_IDX = 0
     SAMPLE_FASTQ_IDX = 2
     _wFragmentRecruitment(sampleReads | map { sample -> [sample[SAMPLE_NAME_IDX], sample[SAMPLE_FASTQ_IDX]] }, genomes)
}


workflow _wFragmentRecruitment {
   take: 
     sampleReads
     genomes
   main:
     fragmentRecruitmentGenomes = params.tempdir + "/fragmentRecruitmentGenomes"
     file(fragmentRecruitmentGenomes).mkdirs()

     genomes | branch { compressed: file(it).name.endsWith(".gz"); uncompressed: !file(it).name.endsWith(".gz") } | set { fileState }
     
     fileState.compressed | pUnzip | set{unzippedGenomes}
     fileState.uncompressed | mix(unzippedGenomes) | collectFile(newLine: true){genome -> ["genomes_list", file(genome).path]} | set { genomesList } 

     fileState.uncompressed | mix(unzippedGenomes) \
       | collectFile(tempDir: fragmentRecruitmentGenomes, sort: params?.steps?.fragmentRecruitment?.frhit?.sort){ genome -> ["genomes",genome.text] } \
       | combine(genomesList) | set { genomesCombined }

     SAMPLE_NAME_IDX = 0
     SAMPLE_FASTQ_IDX = 2
     sampleReads | combine(genomesCombined) | pFrHit

     BAM_FILE_IDX = 1
     pFrHit.out.alignment | map{ align -> align[BAM_FILE_IDX] } | collect | set { alignments }
     pCombinedAlignmentAnalysis(alignments, genomesCombined)
}
