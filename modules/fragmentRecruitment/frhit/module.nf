nextflow.enable.dsl=2
include { wUnmappedReadsFile } from '../../sampleAnalysis/module'

def getOutput(SAMPLE, RUNID, TOOL, filename){
    return SAMPLE + '/' + RUNID + '/' + params.modules.fragmentRecruitment.name + '/' + 
         params.modules.fragmentRecruitment.version.major + "."
         params.modules.fragmentRecruitment.version.minor + "."
         params.modules.fragmentRecruitment.version.patch
          + '/' + TOOL + '/' + filename
}


def getAggregatedOutput(RUNID, TOOL, filename){
    return AGGREGATED + '/' + RUNID + '/' + params.modules.fragmentRecruitment.name + '/' + 
         params.modules.fragmentRecruitment.version.major + "."
         params.modules.fragmentRecruitment.version.minor + "."
         params.modules.fragmentRecruitment.version.patch
          + '/' + TOOL + '/' + filename
}

process pFrHit {

    label 'large'

    tag "$sample"

    stageInMode 'copy'

    publishDir params.output, mode: "${params.publishDirMode}", saveAs: { filename -> getOutput("${sample}", params.runid, "frhit",  filename) }

    container "${params.samtools_bwa_image}"

    when params.steps.containsKey("fragmentRecruitment")

    input:
    tuple val(sample), file(reads), file(genomesCombined), file(genomesList) 

    output:
    tuple val("${sample}"), file("${sample}.bam"), optional: true, emit: alignment
    tuple val("${sample}"), file("coverage"), optional: true, emit: coverageStats
    tuple file(".command.sh"), file(".command.out"), file(".command.err"), file(".command.log")

    shell:
    '''
    mkdir coverage
    seqkit fq2fa !{reads} -o reads.fasta.gz
    AVG_LEN=$(seqkit stats -T reads.fasta.gz | tail -n 1 | cut -f 7)
    HALF_AVG_LEN=$(bc <<< "${AVG_LEN} / 2")
    gunzip reads.fasta.gz
    echo "Min. Coverage Parameter: ${HALF_AVG_LEN}"
    fr-hit !{params.steps.fragmentRecruitment.frhit.additionalParams.frhit} -f 1 -m ${HALF_AVG_LEN} -T !{task.cpus} -a reads.fasta -d !{genomesCombined} -o coverage/out.psl
    if [ -s "coverage/out.psl" ] 
    then
      psl2sam.pl coverage/out.psl | samtools view -bT !{genomesCombined} - | samtools calmd -E - genomes | samtools view -Sb - | samtools sort -o - out > !{sample}.bam 2> /dev/null
      coverm genome -t !{task.cpus} !{params.steps.fragmentRecruitment.frhit.additionalParams.coverm} -b !{sample}.bam --genome-fasta-list !{genomesList} --methods count --output-file coverage/readCount.tsv
      coverm genome -t !{task.cpus} !{params.steps.fragmentRecruitment.frhit.additionalParams.coverm} -b !{sample}.bam --genome-fasta-list !{genomesList} --methods covered_fraction --output-file coverage/coveredFraction.tsv
      coverm genome -t !{task.cpus} !{params.steps.fragmentRecruitment.frhit.additionalParams.coverm} -b !{sample}.bam --genome-fasta-list !{genomesList} --methods covered_bases --output-file coverage/coveredBases.tsv
      coverm genome -t !{task.cpus} !{params.steps.fragmentRecruitment.frhit.additionalParams.coverm} -b !{sample}.bam --genome-fasta-list !{genomesList} --methods length --output-file coverage/genomeLength.tsv
      coverm genome -t !{task.cpus} !{params.steps.fragmentRecruitment.frhit.additionalParams.coverm} -b !{sample}.bam --genome-fasta-list !{genomesList} --methods trimmed_mean --output-file coverage/trimmedMean.tsv
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
  < $x zcat --force > ${x.baseName}
  """
}



process pCombinedAlignmentAnalysis {

    label 'medium'

    when params.steps.containsKey("fragmentRecruitment")

    publishDir params.output, mode: "${params.publishDirMode}", saveAs: { filename -> getAggregatedOutput(params.runid, "frhit",  filename) }

    container "${params.samtools_bwa_image}"

    input:
    file(alignments)
    tuple file(genomesCombined), file(genomesList)

    output:
    path("coverage"), emit: coverageStats
    tuple file(".command.sh"), file(".command.out"), file(".command.err"), file(".command.log")

    shell:
    '''
    mkdir coverage
    samtools merge -@ !{task.cpus} combined_alignments.bam !{alignments}
    coverm genome -t !{task.cpus} !{params.steps.fragmentRecruitment.frhit.additionalParams.coverm} -b combined_alignments.bam --genome-fasta-list !{genomesList} --methods count --output-file coverage/readCount.tsv
    coverm genome -t !{task.cpus} !{params.steps.fragmentRecruitment.frhit.additionalParams.coverm} -b combined_alignments.bam --genome-fasta-list !{genomesList} --methods covered_fraction --output-file coverage/coveredFraction.tsv
    coverm genome -t !{task.cpus} !{params.steps.fragmentRecruitment.frhit.additionalParams.coverm} -b combined_alignments.bam --genome-fasta-list !{genomesList} --methods covered_bases --output-file coverage/coveredBases.tsv
    coverm genome -t !{task.cpus} !{params.steps.fragmentRecruitment.frhit.additionalParams.coverm} -b combined_alignments.bam --genome-fasta-list !{genomesList} --methods length --output-file coverage/genomeLength.tsv
    coverm genome -t !{task.cpus} !{params.steps.fragmentRecruitment.frhit.additionalParams.coverm} -b combined_alignments.bam --genome-fasta-list !{genomesList} --methods trimmed_mean --output-file coverage/trimmedMean.tsv
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

     genomes | splitCsv(sep: '\t', header: true) \
       | map { genome -> file(genome.BINS) } | set {magsList}

     _wFragmentRecruitment(sampleReadsList, magsList)
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

     genomes | pUnzip | set {unzippedGenomes}

     unzippedGenomes | collectFile(newLine: true){genome -> ["genomes_list", file(genome).path]} \
       | set { genomesList } 

     unzippedGenomes \
       | collectFile(tempDir: fragmentRecruitmentGenomes, sort: params?.steps?.fragmentRecruitment?.frhit?.sort){ genome -> ["genomes", genome.text] } \
       | combine(genomesList) | set { genomesCombined }

     SAMPLE_NAME_IDX = 0
     SAMPLE_FASTQ_IDX = 2
     sampleReads | combine(genomesCombined) | pFrHit

     BAM_FILE_IDX = 1
     pFrHit.out.alignment | map{ align -> align[BAM_FILE_IDX] } | collect | set { alignments }
     pCombinedAlignmentAnalysis(alignments, genomesCombined)
}
