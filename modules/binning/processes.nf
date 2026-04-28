
process pGetBinStatistics {

    container "${params.samtools_image}"

    tag "Sample: $sample"

    publishDir params.output, mode: "${params.publishDirMode}", saveAs: { filename -> Output.getOutput("${sample}", params.runid, "${binner}", module, filename) }

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

    publishDir params.output, mode: "${params.publishDirMode}", saveAs: { filename -> Output.getOutput("${sample}", params.runid, "${outputToolDir}", module, filename) }, \
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
	saveAs: { filename -> Output.getOutput("${prefixOutput}", params.runid, "${outputToolDir}", module, filename) }, \
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
    output = Output.getOutput(prefixOutput, params.runid, outputToolDir, module, "")
    template('coverm.sh')
}


process pMinimap2 {

    container "${params.samtools_image}"

    label 'highmemLarge'

    tag "Sample: $sample"

    cache 'deep'

    publishDir params.output, mode: "${params.publishDirMode}", \
	saveAs: { filename -> Output.getOutput("${sample}", params.runid, "${outputToolDir}", module, filename) }

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


process pGetMappingQuality {

    container "${params.samtools_image}"

    tag "Sample: ${sample}"

    publishDir params.output, mode: "${params.publishDirMode}", saveAs: { filename -> Output.getOutput("${sample}", params.runid, "readMappingQuality", module, filename) }

    label 'tiny'

    input:
    val(module)
    tuple val(sample), path(bam)

    output:
    tuple val("${sample}"), file("${sample}_flagstat.tsv"), emit: flagstatRaw
    tuple val("${sample}"), file("${sample}_flagstat_passed.tsv"), emit: flagstatPassed
    tuple val("${sample}"), file("${sample}_flagstat_failed.tsv"), emit: flagstatFailed
    tuple file(".command.sh"), file(".command.out"), file(".command.err"), file(".command.log")

    shell:
    template('mapping_quality.sh')
}


workflow wGetMappingQuality {
    take:
        module
        bam
    main:
        pGetMappingQuality(module, bam)        
    emit:
        flagstatRaw = pGetMappingQuality.out.flagstatRaw
        flagstatPassed = pGetMappingQuality.out.flagstatPassed
        flagstatFailed = pGetMappingQuality.out.flagstatFailed
}

process pBowtie2 {

    container "${params.bowtie_image}"

    label 'highmemLarge'

    tag "Sample: $sample"

    cache 'deep'

    publishDir params.output, mode: "${params.publishDirMode}", \
	saveAs: { filename -> Output.getOutput("${sample}", params.runid, "${outputToolDir}", module, filename) }

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
	saveAs: { filename -> Output.getOutput("${sample}", params.runid, "${outputToolDir}", module, filename) }

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
	saveAs: { filename -> Output.getOutput("${sample}", params.runid, "${outputToolDir}", module, filename) }

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

    memory { Utils.getMemoryResources(params.resources.small, "${sample}", task.attempt, params.resources) }

    cpus { Utils.getCPUsResources(params.resources.small, "${sample}", task.attempt, params.resources) }

    publishDir params.output, mode: "${params.publishDirMode}", \
	saveAs: { filename -> Output.getOutput("${sample}", params.runid, "${outputToolDir}", module, filename) }

    when:
    run

    input:
    val(run)
    tuple val(module), val(outputToolDir), val(metabatParams)
    tuple val(sample), path(contigs), path(bam), val(medianQuality)

    output:
    tuple val("${sample}"), path("${sample}_bin.*.fa", arity: '1..*'), optional: true, emit: bins
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

process pSemiBin2 {

    container "${params.semibin2_image}"

    tag "$sample"

    label 'small'

    memory { Utils.getMemoryResources(params.resources.small, "${sample}", task.attempt, params.resources) }

    cpus { Utils.getCPUsResources(params.resources.small, "${sample}", task.attempt, params.resources) }

    publishDir params.output, mode: "${params.publishDirMode}", \
	saveAs: { filename -> Output.getOutput("${sample}", params.runid, "${outputToolDir}", module, filename) }

    when:
    run

    input:
    val(run)
    tuple val(module), val(outputToolDir), val(semibinParams)
    tuple val(sample), path(contigs), path(bam)

    output:
    tuple val("${sample}"), path("${sample}_bin.*.fa", arity: '1..*'), optional: true, emit: bins
    tuple val("${sample}"), file("${sample}_notBinned.fa"), optional: true, emit: notBinned
    tuple val("${sample}"), file("${sample}_bin_contig_mapping.tsv"), optional: true, emit: binContigMapping
    tuple file(".command.sh"), file(".command.out"), file(".command.err"), file(".command.log")

    script:
    template 'semibin2.sh'
}

/*
* Training for all samples is done seperately. However it uses the input features that account for all the sample data.
*/ 
process pSemiBin2Training {

    container "${params.semibin2_image}"

    tag "Group: $group, Sample: $sample"

    publishDir params.output, mode: "${params.publishDirMode}", \
	saveAs: { filename -> Output.getOutput("${sample}", params.runid, "${outputToolDir}", module, filename) }

    memory { Utils.getMemoryResources(params.resources.medium, "${group}", task.attempt, params.resources) }

    cpus { Utils.getCPUsResources(params.resources.medium, "${group}", task.attempt, params.resources) }

    input:
    tuple val(module), val(outputToolDir), val(semibinParams)
    tuple val(group), val(sample), path(featureOutput), path(featureCompressedOutput)  

    output:
    tuple val("${group}"), val("${sample}"), path("${featureCompressedOutput}", includeInputs: true), optional: true, emit: featureGeneration 
    tuple val("${group}"), val("${sample}"), path("${sample}_output"), optional: true, emit: trainingOutput 
    tuple path(".command.sh"), path(".command.out"), path(".command.err"), path(".command.log")

    script:
    """
    SemiBin2 train_self ${semibinParams} \
        --processes ${task.cpus} \
        --data output/samples/${sample}_contigs/data.csv \
        --data-split output/samples/${sample}_contigs/data_split.csv \
        --output ${sample}_output
    """
}

/*
*
* This process executes the SemiBin2 binning step. 
*
*/
process pSemiBin2Binning {

    container "${params.semibin2_image}"

    tag "Group: $group, Sample: $sample"

    label 'small'

    publishDir params.output, mode: "${params.publishDirMode}", \
	saveAs: { filename -> Output.getOutput("${sample}", params.runid, "${outputToolDir}", module, filename) }

    memory { Utils.getMemoryResources(params.resources.medium, "${group}", task.attempt, params.resources) }

    cpus { Utils.getCPUsResources(params.resources.medium, "${group}", task.attempt, params.resources) }

    input:
    tuple val(module), val(outputToolDir), val(semibinParams)
    val(sequenceType)
    tuple val(group), val(sample), path(trainingOutput), path(contigs), \
        path(featureOutput), path(compressedFeatureOutput), path(sampleGroupsFile) 

    output:
    tuple val("${sample}"), path("${sample}_bin.*.fa", arity: '1..*'), optional: true, emit: bins
    tuple val("${sample}"), file("${sample}_notBinned.fa"), optional: true, emit: notBinned
    tuple val("${sample}"), path("${sampleGroupsFile}", includeInputs: true), optional: true, emit: groupsFile 
    tuple val("${sample}"), file("${sample}_bin_contig_mapping.tsv"), optional: true, emit: binContigMapping
    tuple file(".command.sh"), file(".command.out"), file(".command.err"), file(".command.log")

    script:
    semibin2Command =  sequenceType == "LONG" ? "bin_long" : "bin_short"
    """
    # Run SemiBin2  
    SemiBin2 ${semibin2Command} ${semibinParams} \
        -i ${contigs} \
        --processes ${task.cpus} \
        --compression none \
        --model ${trainingOutput}/model.pt \
        --data ${featureOutput}/samples/${sample}_contigs/data.csv \
        -o binningOutput

    # Process all produced bins 
    BIN_CONTIG_MAPPING=${sample}_bin_contig_mapping.tsv
    echo -e "BIN_ID\tCONTIG\tBINNER" > \${BIN_CONTIG_MAPPING}
    for bin in \$(find binningOutput/output_bins/ -name "*.fa"); do 

	    # Get id of the bin (e.g get 2 of the bin SemiBin_2.fa)
	    ID=\$(echo \$(basename \$bin) | rev | cut -d '_' -f 1 | rev | cut -d '.' -f 1  )
	    BIN_NAME="${sample}_bin.\${ID}.fa"
    
        # Append bin id to every header
	    seqkit replace  -p '(.*)' -r "\\\${1} MAG=\${ID}" \$bin > \${BIN_NAME}

	    # Create bin to contig mapping and add the used binner to each line
	    grep ">" \${bin} | sed 's/>//g' \
		    | sed "s/^/\${BIN_NAME}\t/g;s/\$/\tsemibin2/" >> \${BIN_CONTIG_MAPPING}
    done

    # return not binned fasta files
    BINNED_IDS=binned.tsv
    NOT_BINNED=${sample}_notBinned.fa
    grep -h ">" \$(find binningOutput/output_bins/ -name "*.fa") | tr -d ">" > \${BINNED_IDS}
    if [ -s \${BINNED_IDS} ]; then
	    # Get all not binned Ids
	    seqkit grep -vf \${BINNED_IDS} ${contigs} \
		    | seqkit replace  -p '(.*)' -r "\\\${1} MAG=NotBinned" > \${NOT_BINNED}
    else
	    seqkit replace  -p '(.*)' -r "\\\${1} MAG=NotBinned" ${contigs} > \${NOT_BINNED}
    fi

    # Fix for ownership issue https://github.com/nextflow-io/nextflow/issues/4565
    chmod a+rw -R output binningOutput
    """
}

/*
* This process runs SemiBin2 to collect all sequence related features on the concatenated assembly.
*/
process pSemiBin2GenerateSequenceFeatures {

    container "${params.semibin2_image}"

    tag "Group: $group"

    memory { Utils.getMemoryResources(params.resources.small, "${group}", task.attempt, params.resources) }

    cpus { Utils.getCPUsResources(params.resources.small, "${group}", task.attempt, params.resources) }

    input:
    val(semibinParams)
    tuple val(group), path(bam), path(contigs) 

    output:
    tuple val("${group}"), file("output"), file("feature_generation.tar.gz"), optional: true, emit: features

    script:
    """
    SemiBin2 generate_sequence_features_multi \
    ${semibinParams} \
    --processes ${task.cpus} \
    -i ${contigs} \
    -b ${bam} \
    -o output  

    tar czvf feature_generation.tar.gz .command.sh .command.out .command.err .command.log output
    """
}


/*
* This process saves the information of which sample belongs to which binning group in a file. 
*/
process pSemiBin2SaveGroupingInfo {

    container "${params.semibin2_image}"

    tag "Group: $group"

    label "tiny"

    input:
    tuple val(group), val(samples) 

    output:
    tuple val("${group}"), path("samples_in_group_${group}.tsv") 

    script:
    """
    OUTPUT=samples_in_group_${group}.tsv
    echo -e "GROUP\tSAMPLE" > \${OUTPUT}
    for sample in ${samples.join(' ')}; do
        echo -e "${group}\t\$sample" >> \${OUTPUT}
    done
    """
}

/*
*
* This workflow runs the mapping tool that was selected by the user.
* module: This parameter selects the output path defined by the module
* condition: This parameter contains the tool name and executes the corresponding mapper. 
* option: This parameter defined all parameters that are needed for the tool itself and samtools.
* mapperInput: The contigs and reads that are used as input.
*
*/ 
workflow _wRunMappers {
  take:
    module
    condition
    options
    mapperInput
  main:
    pBowtie2(
        condition | map { mapper -> mapper == "bowtie" },
        options,
        mapperInput
    )

    pBwa(
        condition | map { mapper -> mapper == "bwa" },
        options,
        mapperInput
    )

    pBwa2(
        condition | map { mapper -> mapper == "bwa2" },
        options, 
        mapperInput
    )

    pBowtie2.out.mappedReads 
    | mix(pBwa.out.mappedReads, pBwa2.out.mappedReads) | set { mappedReads }
    pBowtie2.out.unmappedReads 
    | mix(pBwa.out.unmappedReads, pBwa2.out.unmappedReads) | set { unmappedReads }

    wGetMappingQuality(module, mappedReads)
   emit:
     mappedReads = mappedReads
     unmappedReads = unmappedReads
}

/*
* This workflow executes the SemiBin2 binning workflow where the training and binning
* procedure can be parallelized.
*
*/
workflow _wMultiBinningSemiBin2 {
    take:
       sequenceType
       featureGenerationConfig
       trainingConfig
       binningConfig
       sampleGroups
       concatContigs
       sampleContigs
       mapping
    main:
      SAMPLE_IDX = 0
      GROUP_IDX = 0

      // First we try to detect sequence related features on the concatenated contigs.
      sampleGroups
        | combine(mapping, by: SAMPLE_IDX)
        | map { _sample, group, groupSize, bam -> tuple( groupKey(group, groupSize), bam ) }
        | groupTuple(remainder: true)
        | join(concatContigs)
        | set { semibin2FeatureGenerationInput }

     pSemiBin2GenerateSequenceFeatures(featureGenerationConfig, semibin2FeatureGenerationInput)

      // Then we train Semibin on each sample independently using all features detected previously.
     sampleGroups 
       | map{ sample, group, _groupSize -> [group, sample]} 
       | combine(pSemiBin2GenerateSequenceFeatures.out.features, by:GROUP_IDX)
       | set { semibin2TrainingInput } 

     pSemiBin2Training(trainingConfig, semibin2TrainingInput)


     sampleGroups 
        | map { sample, group, groupSize -> tuple( groupKey(group, groupSize), sample ) }
        | groupTuple(remainder: true)
        | pSemiBin2SaveGroupingInfo
        | set { sampleGroupsFile }

     pSemiBin2Training.out.trainingOutput  
     | map { group, sample, trainingOutput -> [sample, group, trainingOutput] }
     | join(sampleContigs, by: SAMPLE_IDX)
     | map { sample, group, trainingOutput, contigs -> [group, sample, trainingOutput, contigs]} 
     | combine(pSemiBin2GenerateSequenceFeatures.out.features, by: GROUP_IDX) 
     | combine(sampleGroupsFile, by: GROUP_IDX) 
     | set { semibin2Input }

    // Now we execute the actual binning step.
     pSemiBin2Binning(binningConfig, sequenceType, semibin2Input)

    emit:
     bins = pSemiBin2Binning.out.bins 
     notBinned = pSemiBin2Binning.out.notBinned 
     binContigMapping = pSemiBin2Binning.out.binContigMapping

}
