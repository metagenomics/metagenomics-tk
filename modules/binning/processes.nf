

String getOutput(SAMPLE, RUNID, MODULE , TOOL, filename){

    def module = ""
    if(MODULE.isEmpty()){
       module = params.modules.binning.name + '/' +
          params.modules.binning.version.major + "." +
          params.modules.binning.version.minor + "." +
          params.modules.binning.version.patch;
    } else {
       module = MODULE
    }

    return SAMPLE + '/' + RUNID + '/' + module  +
          '/' + TOOL + '/' + filename
}

process pGetBinStatistics {

    container "${params.samtools_image}"

    tag "$sample"

    publishDir params.output, mode: "${params.publishDirMode}", saveAs: { filename -> getOutput("${sample}", params.runid, "", "${binner}", filename) }

    label 'tiny'

    input:
    tuple val(sample), path(binContigMapping), path(bam), val(binner), path(bins)

    output:
    tuple val("${sample}"), file("${sample}_contigs_depth.tsv"), optional: true, emit: contigsDepth
    tuple val("${sample}"), file("${sample}_bins_stats.tsv"), optional: true, emit: binsStats
    tuple file(".command.sh"), file(".command.out"), file(".command.err"), file(".command.log")

    shell:
    template 'binStats.sh'
}


process pCovermContigsCoverage {

    label 'medium'

    tag "$sample"

    publishDir params.output, mode: "${params.publishDirMode}", saveAs: { filename -> getOutput("${sample}", params.runid, "${module}" , "${outputToolDir}", filename) }, \
        pattern: "{**.tsv,**.fasta.gz}"


    input:
    val(run)
    tuple val(module), val(outputToolDir), val(covermParams)
    tuple val(sample), path(bamFile)

    when:
    run

    output:
    tuple val("${sample}"), path("${sample}_default_coverm_coverage.tsv"), path("${sample}_metabat_coverm_coverage.tsv"), emit: coverage, optional: true
    tuple val("${sample}"), path("${sample}_default_coverm_coverage.tsv"), emit: default_coverage, optional: true
    tuple val("${sample}"), path("${sample}_metabat_coverm_coverage.tsv"), emit: metabat_coverage, optional: true
    tuple file(".command.out"), file(".command.err"), file(".command.log"), emit: logs

    shell:
    '''
    coverm contig --threads !{task.cpus} \
	--bam-files !{bamFile} !{covermParams} \
	--methods mean trimmed_mean variance length count reads_per_base rpkm tpm \
	| sed '1 s/^/SAMPLE\t/' \
	| sed "2,$ s/^/!{sample}\t/g" > !{sample}_default_coverm_coverage.tsv

    coverm contig --threads !{task.cpus} \
	--bam-files !{bamFile} !{covermParams} \
        --methods metabat \
	| sed '1 s/^/SAMPLE\t/' \
	| sed "2,$ s/^/!{sample}\t/g" > !{sample}_metabat_coverm_coverage.tsv
    '''
}


process pBowtie2 {

    container "${params.bowtie_image}"

    label 'large'

    tag "Sample: $sample"

    publishDir params.output, mode: "${params.publishDirMode}", saveAs: { filename -> getOutput("${sample}", params.runid, "${module}", "${outputToolDir}", filename) }

    input:
    val(run)
    tuple val(module), val(outputToolDir), val(bowtieParams), val(getUnmapped)
    tuple val(sample), path(contigs), path(pairedReads, stageAs: 'paired.fq.gz'), path(unpairedReads, stageAs: 'unpaired.fq.gz')

    when:
    run

    output:
    tuple val("${sample}"), file("${sample}.bam"), optional: true, emit: mappedReads
    tuple val("${sample}"), file("${sample}_unmapped.fq.gz"), optional: true, emit: unmappedReads
    tuple val("${sample}"), file("${sample}_bowtie_stats.txt"), optional: true, emit: stats
    tuple file(".command.sh"), file(".command.out"), file(".command.err"), file(".command.log")

    shell:
    getUnmapped = getUnmapped ? "TRUE" : ""
    '''
    INDEX=!{sample}.index
    # Build Bowtie Index
    bowtie2-build --threads !{task.cpus} --quiet !{contigs} $INDEX

    # Run Bowtie
    bowtie2 -p !{task.cpus} !{bowtieParams} -x $INDEX \
              --interleaved paired.fq.gz -U unpaired.fq.gz 2> !{sample}_bowtie_stats.txt \
             | samtools view -F 3584 --threads !{task.cpus} -bS - \
             | samtools sort -l 9 --threads !{task.cpus} - > !{sample}.bam

    # If Fragment Recruitment is selected then reads that could not be mapped should be returned
    if [[ "!{getUnmapped}" == "TRUE" ]]; then
        samtools bam2fq -f 4 !{sample}.bam | pigz --best --processes !{task.cpus} > !{sample}_unmapped.fq.gz
    fi
    '''
}

