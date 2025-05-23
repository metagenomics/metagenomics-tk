

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

    tag "Sample: $sample"

    publishDir params.output, mode: "${params.publishDirMode}", saveAs: { filename -> getOutput("${sample}", params.runid, "${module}", "${binner}", filename) }

    memory { Utils.getMemoryResources(params.resources.tiny, "${sample}", task.attempt, params.resources) }

    cpus { Utils.getCPUsResources(params.resources.tiny, "${sample}", task.attempt, params.resources) }

    input:
    val(module)
    tuple val(sample), path(binContigMapping), path(bam), val(binner), path(bins), val(medianQuality)

    output:
    tuple val("${sample}"), file("${sample}_contigs_depth.tsv"), optional: true, emit: contigsDepth
    tuple val("${sample}"), file("${sample}_bins_stats.tsv"), optional: true, emit: binsStats
    tuple file(".command.sh"), file(".command.out"), file(".command.err"), file(".command.log")

    shell:
    DO_NOT_ESTIMATE_QUALITY = -1 
    MEDIAN_QUALITY=Double.parseDouble(medianQuality)
    percentIdentity = MEDIAN_QUALITY != DO_NOT_ESTIMATE_QUALITY ? \
	" --percentIdentity="+Utils.getMappingIdentityParam(MEDIAN_QUALITY) : " "

    template 'binStats.sh'
}


process pCovermContigsCoverage {

    label 'small'

    tag "Sample: $sample"

    publishDir params.output, mode: "${params.publishDirMode}", saveAs: { filename -> getOutput("${sample}", params.runid, "${module}" , "${outputToolDir}", filename) }, \
        pattern: "{**.tsv,**.fasta.gz,**.out,**.err,**.log,**.sh}"

    input:
    val(run)
    tuple val(module), val(outputToolDir), val(covermParams)
    tuple val(sample), path(bamFile), val(medianQuality)

    when:
    run

    output:
    tuple val("${sample}"), path("${sample}_default_coverm_coverage.tsv"), path("${sample}_metabat_coverm_coverage.tsv"), emit: coverage, optional: true
    tuple val("${sample}"), path("${sample}_default_coverm_coverage.tsv"), emit: default_coverage, optional: true
    tuple val("${sample}"), path("${sample}_metabat_coverm_coverage.tsv"), emit: metabat_coverage, optional: true
    tuple file(".command.out"), file(".command.err"), file(".command.log"), file(".command.sh"), emit: logs

    shell:
    DO_NOT_ESTIMATE_QUALITY = -1 
    MEDIAN_QUALITY=Double.parseDouble(medianQuality)
    percentIdentity = MEDIAN_QUALITY != DO_NOT_ESTIMATE_QUALITY ? \
	" --min-read-percent-identity "+Utils.getMappingIdentityParam(MEDIAN_QUALITY) : " "
    '''
    coverm contig --threads !{task.cpus} \
	--bam-files !{bamFile} !{covermParams} !{percentIdentity} \
	--methods mean trimmed_mean variance length count reads_per_base rpkm tpm \
	| sed '1 s/^.*$/SAMPLE\tCONTIG\tMEAN_CONTIG\tTRIMMED_MEAN_CONTIG\tVARIANCE_CONTIG\tLENGTH_CONTIG\tREAD_COUNT_CONTIG\tREADS_PER_BASE_CONTIG\tRPKM_CONTIG\tTPM_CONTIG/g' \
	| sed "2,$ s/^/!{sample}\t/g" > !{sample}_default_coverm_coverage.tsv

    coverm contig --threads !{task.cpus} \
	--bam-files !{bamFile} !{covermParams} !{percentIdentity} \
        --methods metabat \
	| sed '1 s/^.*$/SAMPLE\tCONTIG_NAME\tMETABAT_CONTIG_LENGTH\tMETABAT_TOTAL_AVG_DEPTH_CONTIG\tMETABAT_BAM_CONTIG\tMETABAT_VARIANCE_CONTIG/g' \
	| sed "2,$ s/^/!{sample}\t/g" > !{sample}_metabat_coverm_coverage.tsv
    '''
}



process pCovermGenomeCoverage {

    label 'small'

    tag "Sample: $sample"

    publishDir params.output, mode: "${params.publishDirMode}", \
	saveAs: { filename -> getOutput("${prefixOutput}", params.runid, "${module}" , "${outputToolDir}", filename) }, \
        pattern: "{**.tsv}"

    input:
    val(run)
    val(prefix)
    tuple val(module), val(outputToolDir), val(covermParams)
    tuple val(sample), file(mapping), file(index), file(list_of_representatives), val(medianQuality)

    when:
    run

    output:
    tuple val("${sample}"), path("${sample}_mean.tsv"), emit: mean
    tuple val("${sample}"), path("${sample}_trimmed_mean.tsv"), emit: trimmedMean
    tuple val("${sample}"), path("${sample}_count.tsv"), emit: count
    tuple val("${sample}"), path("${sample}_rpkm.tsv"), emit: rpkm
    tuple val("${sample}"), path("${sample}_relative_abundance.tsv"), emit: relativeAbundance
    tuple val("${sample}"), path("${sample}_tpm.tsv"), emit: tpm
    tuple val("${sample}"), val(output), \
	val(params.LOG_LEVELS.INFO), file(".command.sh"), \
	file(".command.out"), file(".command.err"), file(".command.log"), emit: logs

    shell:
    DO_NOT_ESTIMATE_QUALITY = -1 
    MEDIAN_QUALITY=Double.parseDouble(medianQuality)
    percentIdentity = MEDIAN_QUALITY != DO_NOT_ESTIMATE_QUALITY ? \
	" --min-read-percent-identity "+Utils.getMappingIdentityParam(MEDIAN_QUALITY) : " "
    prefixOutput = prefix.trim() ? prefix : sample
    output = getOutput(prefixOutput, params.runid, module, outputToolDir, "")
    template('coverm.sh')
}


process pMinimap2 {

    container "${params.samtools_image}"

    label 'highmemLarge'

    tag "Sample: $sample"

    cache 'deep'

    publishDir params.output, mode: "${params.publishDirMode}", \
	saveAs: { filename -> getOutput("${sample}", params.runid, "${module}", "${outputToolDir}", filename) }

    input:
    val(run)
    tuple val(module), val(outputToolDir), val(minimapParams), val(samtoolsViewParams), val(getUnmapped)
    tuple val(sample), path(contigs), path(reads, stageAs: 'reads.fq.gz')

    when:
    run

    output:
    tuple val("${sample}"), file("${sample}.bam"), optional: true, emit: mappedReads
    tuple val("${sample}"), file("${sample}_unmapped.fq.gz"), optional: true, emit: unmappedReads
    tuple file(".command.sh"), file(".command.out"), file(".command.err"), file(".command.log")

    shell:
    getUnmapped = getUnmapped ? "TRUE" : ""
    '''
    minimap2 -t !{task.cpus} !{minimapParams} -ax map-ont !{contigs} reads.fq.gz \
             | samtools view !{samtoolsViewParams} --threads !{task.cpus} -bS - \
             | samtools sort -l 9 --threads !{task.cpus} - > !{sample}.bam

    # If Fragment Recruitment is selected then reads that could not be mapped should be returned
    if [[ "!{getUnmapped}" == "TRUE" ]]; then
        samtools bam2fq -f 4 !{sample}.bam | pigz --best --processes !{task.cpus} > !{sample}_unmapped.fq.gz
    fi
    '''
}



process pBowtie2 {

    container "${params.bowtie_image}"

    label 'highmemLarge'

    tag "Sample: $sample"

    cache 'deep'

    publishDir params.output, mode: "${params.publishDirMode}", \
	saveAs: { filename -> getOutput("${sample}", params.runid, "${module}", "${outputToolDir}", filename) }

    input:
    val(run)
    tuple val(module), val(outputToolDir), val(bowtieParams), val(samtoolsViewParams), val(getUnmapped)
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
             | samtools view !{samtoolsViewParams} --threads !{task.cpus} -bS - \
             | samtools sort -l 9 --threads !{task.cpus} - > !{sample}.bam

    # If Fragment Recruitment is selected then reads that could not be mapped should be returned
    if [[ "!{getUnmapped}" == "TRUE" ]]; then
        samtools bam2fq -f 4 !{sample}.bam | pigz --best --processes !{task.cpus} > !{sample}_unmapped.fq.gz
    fi
    '''
}


process pBwa {

    container "${params.samtools_bwa_image}"

    label 'highmemLarge'

    cache 'deep'

    tag "Sample: $sample"

    publishDir params.output, mode: "${params.publishDirMode}", \
	saveAs: { filename -> getOutput("${sample}", params.runid, "${module}", "${outputToolDir}", filename) }

    input:
    val(run)
    tuple val(module), val(outputToolDir), val(bwaParams), val(samtoolsViewParams), val(getUnmapped)
    tuple val(sample), path(contigs), path(pairedReads, stageAs: 'paired.fq.gz'), path(unpairedReads, stageAs: 'unpaired.fq.gz')

    when:
    run

    output:
    tuple val("${sample}"), file("${sample}.bam"), optional: true, emit: mappedReads
    tuple val("${sample}"), file("${sample}_unmapped.fq.gz"), optional: true, emit: unmappedReads
    tuple file(".command.sh"), file(".command.out"), file(".command.err"), file(".command.log")

    shell:
    getUnmapped = getUnmapped ? "TRUE" : ""
    '''
    # Build BWA Index
    bwa index !{contigs}

    # Run BWA
    bwa mem !{bwaParams} -p  \
       -t !{task.cpus} !{contigs} <(cat !{pairedReads} !{unpairedReads}) - \
      | samtools view !{samtoolsViewParams} -@ !{task.cpus} -S -b - \
      | samtools sort -l 9 -@ !{task.cpus} - > !{sample}.bam

    # If Fragment Recruitment is selected then reads that could not be mapped should be returned
    if [[ "!{getUnmapped}" == "TRUE" ]]; then
        samtools bam2fq -f 4 !{sample}.bam | pigz --best --processes !{task.cpus} > !{sample}_unmapped.fq.gz
    fi
    '''
}


process pBwa2 {

    container "${params.samtools_bwa2_image}"

    label 'highmemLarge'

    cache 'deep'

    tag "Sample: $sample"

    publishDir params.output, mode: "${params.publishDirMode}", \
	saveAs: { filename -> getOutput("${sample}", params.runid, "${module}", "${outputToolDir}", filename) }

    input:
    val(run)
    tuple val(module), val(outputToolDir), val(bwaParams), val(samtoolsViewParams), val(getUnmapped)
    tuple val(sample), path(contigs), path(pairedReads, stageAs: 'paired.fq.gz'), path(unpairedReads, stageAs: 'unpaired.fq.gz')

    when:
    run

    output:
    tuple val("${sample}"), file("${sample}.bam"), optional: true, emit: mappedReads
    tuple val("${sample}"), file("${sample}_unmapped.fq.gz"), optional: true, emit: unmappedReads
    tuple file(".command.sh"), file(".command.out"), file(".command.err"), file(".command.log")

    shell:
    getUnmapped = getUnmapped ? "TRUE" : ""
    '''
    # Build BWA2 Index
    bwa-mem2 index !{contigs}

    # Run BWA
    bwa-mem2 mem !{bwaParams} -p  \
       -t !{task.cpus} !{contigs} <(cat !{pairedReads} !{unpairedReads}) - \
      | samtools view !{samtoolsViewParams} -@ !{task.cpus} -S -b - \
      | samtools sort -l 9 -@ !{task.cpus} - > !{sample}.bam

    # If Fragment Recruitment is selected then reads that could not be mapped should be returned
    if [[ "!{getUnmapped}" == "TRUE" ]]; then
        samtools bam2fq -f 4 !{sample}.bam | pigz --best --processes !{task.cpus} > !{sample}_unmapped.fq.gz
    fi
    '''
}


process pMetabat {

    container "${params.metabat_image}"

    tag "$sample"

    label 'small'

    publishDir params.output, mode: "${params.publishDirMode}", \
	saveAs: { filename -> getOutput("${sample}", params.runid, "${module}", "${outputToolDir}", filename) }

    when:
    run

    input:
    val(run)
    tuple val(module), val(outputToolDir), val(metabatParams)
    tuple val(sample), path(contigs), path(bam), val(medianQuality)

    output:
    tuple val("${sample}"), file("${sample}_bin.*.fa"), optional: true, emit: bins
    tuple val("${sample}"), file("${sample}_notBinned.fa"), optional: true, emit: notBinned
    tuple val("${sample}"), file("${sample}_bin_contig_mapping.tsv"), optional: true, emit: binContigMapping
    tuple file(".command.sh"), file(".command.out"), file(".command.err"), file(".command.log")


    shell:
    DO_NOT_ESTIMATE_QUALITY = -1 
    MEDIAN_QUALITY=Double.parseDouble(medianQuality)
    percentIdentity = MEDIAN_QUALITY != DO_NOT_ESTIMATE_QUALITY ? \
	" PCTID="+Utils.getMappingIdentityParam(MEDIAN_QUALITY) : " "
    template 'metabat.sh'
}


