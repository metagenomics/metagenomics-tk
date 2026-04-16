include {
    pGetBinStatistics as pGetBinStatistics ;
    pCovermContigsCoverage ;
    pCovermGenomeCoverage ;
    pBowtie2 ;
    pSemiBin2 ;
    pMetabat ;
    pMinimap2 ;
    pBwa ;
    wGetMappingQuality ;
    _wMultiMappingSemiBin2 ;
    pBwa2 ;
    _wRunMappers;
} from './processes'
include { pProdigal ; pHmmSearch } from '../annotation/module'
include { wSaveSettingsList } from '../config/module'
include { createMap; mapJoin } from '../utils/methods'


include { pDumpLogs } from '../utils/processes'


/*
*
* Semibin2 expects all contig names to have a prefix with the sample name.
*
*/
process pConcatContigs {

    container "${params.semibin2_image}"

    tag "Group: ${group}"

    label 'medium'

    containerOptions params.apptainer ? "" : ' --user 0:0 '

    input:
    tuple val(group), path(contigs, stageAs: 'input_*'), val(samples)

    output:
    tuple val("${group}"), file("${group}.fa.gz"),  emit: contigs

    script:
    def file_list = contigs.join(' ')
    def name_list = samples.join(' ')
    """
    # Convert space-separated strings into Bash arrays
    raw_files=($file_list)
    names=($name_list)

    mkdir contigs
    # Loop through and create symlinks using the sample name
    for i in \${!raw_files[@]}; do
        REAL_PATH=\$(readlink -f \${raw_files[\$i]})
        ln -s "\$REAL_PATH" contigs/\${names[\$i]}_contigs.fa.gz
    done 

    SemiBin2 concatenate_fasta \
    --input-fasta contigs/* \
    --output output

    mv output/concatenated.fa.gz ${group}.fa.gz
    chmod a+rwx ${group}.fa.gz
    """
}

/*
*
* This method renames contigs in the bam file to their original name before they were renamed by SemiBin2.
*
*/
process pRenameContigsInMapping {

    container "${params.samtools_image}"

    tag "Sample: ${sample}"

    label 'small'

    containerOptions params.apptainer ? "" : ' --user 0:0 '

    input:
    tuple val(sample), path(mapping)

    output:
    tuple val("${sample}"), path("output/output.bam"), emit: mapping

    script:
    """
    mkdir output
    samtools view -H ${mapping} \
        | sed "s/${sample}_contigs://g" > new_header.sam
    samtools reheader new_header.sam ${mapping} > output/output.bam
    """
}

/*
 * This entrypoint takes the following tab separated file of files containing contigs, group information for co binning, paired and optional single reads 
 * as input and produces binning results.
 * The input file for reads must have two columns seperated by tabs:
 * SAMPLE and READS
 * The input file for contigs must have two columns seperated by tabs:
 * SAMPLE and CONTIGS
 * The input file for groups contains two columns: SAMPLE and GROUP.
*/ 
workflow wMultiBinningShortReadFile {
    main:
       SAMPLE_IDX = 0       
       SAMPLE_PAIRED_IDX = 1
       UNPAIRED_IDX = 2
       GROUP_IDX = 1       

       channel.from(file(params.steps.multiBinning.input.contigs)) | splitCsv(sep: '\t', header: true) \
             | map { it -> [ it.SAMPLE, file(it.CONTIGS)]} | set { contigs }

       channel.from(file(params.steps.multiBinning.input.groups)) | splitCsv(sep: '\t', header: true) \
             | map { it -> [ it.SAMPLE, it.GROUP]} | set { groupNames }

       groupNames | groupTuple(by: GROUP_IDX) | map { samplesList, group -> [samplesList.size(), group] }
             | combine(groupNames, by: GROUP_IDX) 
             | map { group, groupSize, sample -> [sample, group, groupSize]} | set { groups }

       readsPaired = channel.empty()
       if(params.steps.multiBinning.input.containsKey("paired")) {
       	 channel.from(file(params.steps.multiBinning.input.paired)) | splitCsv(sep: '\t', header: true) \
             | map { it -> [ it.SAMPLE, file(it.READS)]} | set { readsPaired }
       }

       readsSingle = channel.empty()
       if(params.steps.multiBinning.input.containsKey("single")) {
         channel.from(file(params.steps.multiBinning.input.single)) | splitCsv(sep: '\t', header: true) \
             | map { it -> [ it.SAMPLE, file(it.READS)]} | set { readsSingle }
       }

       readsPaired | join(readsSingle, by: SAMPLE_IDX, remainder: true)
            | map { sample -> sample[UNPAIRED_IDX] == null ? \
		[sample[SAMPLE_IDX], sample[SAMPLE_PAIRED_IDX], file("NOT_SET")] : sample }
		    | set { reads }

       wSaveSettingsList(reads | map { it -> it[SAMPLE_IDX] })

       _wBinningShortRead(contigs, reads, groups)
}

/*
* The format of the channels is described in the corresponding file entrypoint.
*/
workflow wMultiBinningShortReadList {
    take:
    contigs
    inputReads
    sampleGroups

    main:
    _wBinningShortRead(contigs, inputReads, sampleGroups)

    emit:
    binsStats = _wBinningShortRead.out.binsStats
    bins = _wBinningShortRead.out.bins
    mapping = _wBinningShortRead.out.mapping
    notBinnedContigs = _wBinningShortRead.out.notBinnedContigs
    unmappedReads = _wBinningShortRead.out.unmappedReads
    contigCoverage = _wBinningShortRead.out.contigCoverage
}


/*
 * This entrypoint takes the following file of files as input:
 * Contigs: The contigs file must have the columns SAMPLE and CONTIGS where CONTIGS points to the actual file containing contigs.
 * Reads: The reads file must have the columns SAMPLE and READS where READS points to the actual fastq file.
 * Quality: The quality tsv file contains the columns SAMPLE and QUALITY. The qualiy column contains the median quality for every sample.
 * The input file for groups contains two columns: SAMPLE and GROUP.
 */
workflow wMultiBinningLongReadFile {
    main:
       SAMPLE_IDX = 0       
       GROUP_IDX = 1       

       channel.from(file(params.steps.multiBinningONT.input.contigs)) | splitCsv(sep: '\t', header: true) \
             | map { it -> [ it.SAMPLE, file(it.CONTIGS)]} | set { contigs }

       channel.from(file(params.steps.multiBinningONT.input.groups)) | splitCsv(sep: '\t', header: true) \
             | map { it -> [ it.SAMPLE, it.GROUP]} | set { groupNames }

       groupNames | groupTuple(by: GROUP_IDX) | map { samplesList, group -> [samplesList.size(), group] }
             | combine(groupNames, by: GROUP_IDX) 
             | map { group, groupSize, sample -> [sample, group, groupSize]} | set { groups }

       reads = channel.empty()
       if(params.steps.multiBinningONT.input.containsKey("reads")) {
       	 channel.from(file(params.steps.multiBinningONT.input.reads)) | splitCsv(sep: '\t', header: true) \
             | map { it -> [ it.SAMPLE, file(it.READS)]} | set { reads }
       }

       quality = channel.empty()
       if(params.steps.multiBinningONT.input.containsKey("quality")) {
       	 channel.from(file(params.steps.multiBinningONT.input.quality)) | splitCsv(sep: '\t', header: true) \
             | map { it -> [ it.SAMPLE, it.QUALITY]} | set { quality  }
       }

       wSaveSettingsList(reads | map { it -> it[SAMPLE_IDX] })

      _wBinningLongRead(contigs, reads, groups, quality)
}


/*
* The format of the channels is described in the corresponding file entrypoint.
*/
workflow wMultiBinningLongReadList {
    take:
    contigs
    inputReads
    sampleGroups
    medianQuality

    main:
    _wBinningLongRead(contigs, inputReads, sampleGroups, medianQuality)

    emit:
    binsStats = _wBinningLongRead.out.binsStats
    bins = _wBinningLongRead.out.bins
    mapping = _wBinningLongRead.out.mapping
    notBinnedContigs = _wBinningLongRead.out.notBinnedContigs
    unmappedReads = _wBinningLongRead.out.unmappedReads
    contigCoverage = _wBinningLongRead.out.contigCoverage
}

/*
*
* This workflow takes an input_reads channel as input with the following format [SAMPLE, READS PAIRED, READS UNPAIRED]
*
*/
workflow _wBinningLongRead {
    take:
    contigs
    inputReads
    sampleGroups
    medianQuality

    main:
    DO_NOT_ESTIMATE_IDENTITY = "-1"
    SAMPLE_IDX = 0
    GROUP_IDX = 0

    // Concatenate all assemblies per groups using the semibin specific method.
    sampleGroups
        | join(contigs, by: SAMPLE_IDX)
        | map { sample, group, groupSize, filePath ->  tuple( groupKey(group, groupSize), filePath, sample) }
        | groupTuple(remainder: true)
        | set { assemblyList }

    pConcatContigs(assemblyList)

    // Reads of all samples per group must be mapped to the concatenated assembly.
    sampleGroups
        | join(inputReads, by: GROUP_IDX)
        | map { sample, group, _groupSize, readsFilePath -> [group, sample, readsFilePath] }
        | set { readsList }

    pConcatContigs.out.contigs
        | combine(readsList, by: GROUP_IDX)
        | map { _group, concatContigs, sample, pairedFilePath -> [sample, concatContigs, pairedFilePath] }
        | set { mapperInput }

    def modulePath = Output.getModulePath(params.modules.multiBinningONT)

    pMinimap2(
        channel.value(params?.steps?.containsKey("multiBinningONT")),
        channel.value(
            [
                modulePath,
                "contigMapping",
                params.steps?.multiBinningONT?.minimap?.additionalParams?.minimap,
                params.steps?.multiBinningONT?.minimap?.additionalParams?.samtoolsView,
                params.steps.containsKey("fragmentRecruitment"),
            ]
        ),
        mapperInput,
    )

    pMinimap2.out.mappedReads | set { mappedReads }

    mappedReads | wGetMappingQuality

    // We need to rename contigs in the bam file to their original name for consistency in follow up methods.
    mappedReads
        | pRenameContigsInMapping
        | set { renamedContigsMappedReads }

    pCovermContigsCoverage(
        channel.value(params?.steps?.multiBinningONT.find { it.key == "contigsCoverage" }?.value),
        channel.value(
            [
                modulePath,
                "contigCoverage",
                params?.steps?.multiBinningONT?.contigsCoverage?.additionalParams,
            ]
        ),
        renamedContigsMappedReads | join(medianQuality, by: SAMPLE_IDX),
    )

    // Run SembiBin2 workflow 
    _wMultiMappingSemiBin2(sampleGroups, pConcatContigs.out.contigs, contigs, mappedReads)

    _wMultiMappingSemiBin2.out.bins | set { bins }

    _wMultiMappingSemiBin2.out.notBinned | set { notBinned }

    emptyFile = file(params.tempdir + "/empty")
    emptyFile.text = ""

    ALIGNMENT_INDEX = 2
    pCovermGenomeCoverage(
        channel.value(params?.steps?.multiBinningONT.find { it.key == "genomeCoverage" }?.value),
        channel.value(""),
        channel.value(
            [
                modulePath,
                "genomeCoverage",
                params?.steps?.multiBinningONT?.genomeCoverage?.additionalParams,
            ]
        ),
        renamedContigsMappedReads | join(bins, by: SAMPLE_IDX) | map { sample ->
            sample.addAll(ALIGNMENT_INDEX, emptyFile)
            sample
        } | join(medianQuality, by: SAMPLE_IDX),
    )

    pCovermGenomeCoverage.out.logs | pDumpLogs

    // Flatten metabat outputs per sample and create a map with the 
    // following entries [BIN_ID:bin.name, SAMPLE:sample, PATH:bin]
    bins | map { it -> Utils.flattenTuple(it) } | flatMap { it -> createMap(it) } | set { binMap }

    _wMultiMappingSemiBin2.out.binContigMapping
        | join(renamedContigsMappedReads, by: SAMPLE_IDX)
        | combine(channel.from("semibin2"))
        | join(bins, by: SAMPLE_IDX)
        | set { semibin2BinStatisticsInput }

    semibin2BinStatisticsInput
        | combine(channel.value(DO_NOT_ESTIMATE_IDENTITY))
        | set { binStatsInput }

    pGetBinStatistics(channel.value(modulePath), binStatsInput)

    // Add bin statistics 
    pGetBinStatistics.out.binsStats
        | map { _sample, stats -> file(stats) }
        | map { it -> file(it[1]) }
        | splitCsv(sep: '\t', header: true)
        | set { binsStats }
    mapJoin(binsStats, binMap, "BIN_ID", "BIN_ID") | set { binMap }

    emit:
    binsStats = binMap
    bins = bins
    mapping = renamedContigsMappedReads
    notBinnedContigs = notBinned
    unmappedReads =  pMinimap2.out.unmappedReads
    contigCoverage = pCovermContigsCoverage.out.coverage
}


/*
*
* This workflow takes an input_reads channel as input with the following format [SAMPLE, READS PAIRED, READS UNPAIRED]
*
*/
workflow _wBinningShortRead {
    take:
    contigs
    inputReads
    sampleGroups

    main:
    DO_NOT_ESTIMATE_IDENTITY = "-1"
    SAMPLE_IDX = 0
    GROUP_IDX = 0

    // Concatenate all assemblies per groups using the semibin specific method.
    sampleGroups
        | join(contigs, by: SAMPLE_IDX)
        | map { sample, group, groupSize, filePath -> tuple( groupKey(group, groupSize), filePath, sample) }
        | groupTuple(remainder: true)
        | set { assemblyList }

    pConcatContigs(assemblyList)

    // Reads of all samples per group must be mapped to the concatenated assembly.
    sampleGroups
        | join(inputReads, by: GROUP_IDX)
        | map { sample, group, pairedFilePath, unpairedFilePath -> [group, sample, pairedFilePath, unpairedFilePath] }
        | set { readsList }

    pConcatContigs.out.contigs
        | combine(readsList, by: GROUP_IDX)
        | map { _group, concatContigs, sample, pairedFilePath, unpairedFilePath -> [sample, concatContigs, pairedFilePath, unpairedFilePath] }
        | set { mapperInput }


    mappedReads = channel.empty()
    unmappedReads = channel.empty()
    def modulePath = ""
    if (params.steps.containsKey("multiBinning")) {
        // 2. Identify which tool is active (e.g., via a param or checking keys)
        def activeMapper = params.steps?.multiBinning?.keySet()?.find { tool -> ["bwa", "bowtie", "bwa2"].contains(tool) }
        def mapper = params.steps?.multiBinning[activeMapper]
        modulePath = Output.getModulePath(params.modules.multiBinning)

        // 3. Create the dynamic channel
        mapperConfig = mapper
            ? channel.value(
                [
                    modulePath,
                    "contigMapping",
                    mapper.additionalParams[activeMapper],
                    mapper.additionalParams.samtoolsView,
                    params.steps.containsKey("fragmentRecruitment"),
                ]
            )
            : channel.empty()

        _wRunMappers(channel.value(activeMapper), mapperConfig, mapperInput)

        _wRunMappers.out.mappedReads | set { mappedReads }

        _wRunMappers.out.unmappedReads | set { unmappedReads }
    }

    // We need to rename contigs in the bam file to their original name for consistency for follow up methods.
    mappedReads
        | pRenameContigsInMapping
        | set { renamedContigsMappedReads }

    pCovermContigsCoverage(
        channel.value(params?.steps?.multiBinning.find { mapperVal -> mapperVal.key == "contigsCoverage" }?.value),
        channel.value(
            [
                modulePath,
                "contigCoverage",
                params?.steps?.multiBinning?.contigsCoverage?.additionalParams,
            ]
        ),
        renamedContigsMappedReads | combine(channel.value(DO_NOT_ESTIMATE_IDENTITY)),
    )

    // Run SembiBin2 workflow 
    _wMultiMappingSemiBin2(sampleGroups, pConcatContigs.out.contigs, contigs, mappedReads)

    _wMultiMappingSemiBin2.out.bins | set { bins }

    _wMultiMappingSemiBin2.out.notBinned | set { notBinned }

    emptyFile = file(params.tempdir + "/empty")

    ALIGNMENT_INDEX = 2
    pCovermGenomeCoverage(
        channel.value(params?.steps?.multiBinning.find { mapperVal -> mapperVal.key == "genomeCoverage" }?.value),
        channel.value(""),
        channel.value(
            [
                modulePath,
                "genomeCoverage",
                params?.steps?.multiBinning?.genomeCoverage?.additionalParams,
            ]
        ),
        renamedContigsMappedReads | join(bins, by: SAMPLE_IDX) | map { sample ->
            sample.addAll(ALIGNMENT_INDEX, emptyFile)
            sample
        } | combine(channel.value(DO_NOT_ESTIMATE_IDENTITY)),
    )

    pCovermGenomeCoverage.out.logs | pDumpLogs

    // Flatten metabat outputs per sample and create a map with the 
    // following entries [BIN_ID:bin.name, SAMPLE:sample, PATH:bin]
    bins | map { it -> Utils.flattenTuple(it) } | flatMap { it -> createMap(it) } | set { binMap }

    _wMultiMappingSemiBin2.out.binContigMapping
        | join(renamedContigsMappedReads, by: SAMPLE_IDX)
        | combine(channel.from("semibin2"))
        | join(_wMultiMappingSemiBin2.out.bins, by: SAMPLE_IDX)
        | set { semibin2BinStatisticsInput }

    semibin2BinStatisticsInput
        | combine(channel.value(DO_NOT_ESTIMATE_IDENTITY))
        | set { binStatsInput }

    pGetBinStatistics(channel.value(modulePath), binStatsInput)

    // Add bin statistics 
    pGetBinStatistics.out.binsStats
        | map { _sample, stats -> file(stats) }
        | splitCsv(sep: '\t', header: true)
        | set { binsStats }
    mapJoin(binsStats, binMap, "BIN_ID", "BIN_ID") | set { binMap }

    emit:
    binsStats = binMap
    bins = bins
    mapping = renamedContigsMappedReads
    notBinnedContigs = notBinned
    unmappedReads = unmappedReads
    contigCoverage = pCovermContigsCoverage.out.coverage
}
