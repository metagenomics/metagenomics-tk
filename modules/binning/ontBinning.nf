include { wSaveSettingsList } from '../config/module'

include {
    pGetBinStatistics ;
    pCovermContigsCoverage ;
    pCovermGenomeCoverage ;
    pMinimap2 ;
    pMetabat ;
    pSemiBin2 ;
    wGetMappingQuality
} from './processes'

include { createMap ; mapJoin } from '../utils/methods'

include { pDumpLogs } from '../utils/processes'

/*
*
* Binning Metagenomic Contigs via Composition, Coverage and Assembly Graphs. 
* Since MetaCoag is only used for processing Metaflye assemblies, it takes Flyes assemblyInfo as input. 
*
*/
process pMetaCoAG {

    container "${params.metacoag_image}"

    containerOptions params.apptainer ? "" : ' --user 0:0 '

    tag "Sample: ${sample}"

    label 'highmemLarge'

    publishDir params.output, mode: "${params.publishDirMode}", saveAs: { filename -> Output.getOutput("${sample}", params.runid, "metacoag", params.modules.binningONT, filename) }

    when params.steps.containsKey("binningONT") && params.steps.binningONT.containsKey("metacoag")

    input:
    tuple val(sample), path(graph), path(contigs), path(bam), path(headerMapping), path(flyeAssemblyInfo)

    output:
    tuple val("${sample}"), path("${sample}_bin.*.fa", arity: '1..*'), optional: true, emit: bins
    tuple val("${sample}"), file("${sample}_notBinned.fa"), optional: true, emit: notBinned
    tuple val("${sample}"), file("${sample}_bin_contig_mapping.tsv"), optional: true, emit: binContigMapping
    tuple file(".command.sh"), file(".command.out"), file(".command.err"), file(".command.log")

    shell:
    template("metacoag.sh")
}


/*
 * This entrypoint takes the following file of files as input:
 * Contigs: The contigs file must have the columns SAMPLE and CONTIGS where CONTIGS points to the actual file containing contigs.
 * Reads: The reads file must have the columns SAMPLE and READS where READS points to the actual fastq file.
 * Graph: The tsv file must have two columns: SAMPLE and GRAPH. The assembly graph must be in gfa format. 
 * Header Mapping: The tsv file have the columns SAMPLE and MAPPING where the mapping file maps contig names
 * in the assembly to the contig names in the graph file.
 * Assembly Info: The tsv file must have the columns SAMPLE and INFO. Info points to a file that describes the path in the assembly graph
 * for every contig.
 * Quality: The quality tsv file contains the columns SAMPLE and QUALITY. The qualiy column contains the median quality for every sample.
 */
workflow wOntBinningFile {
    SAMPLE_IDX = 0

    channel.from(file(params.steps.binningONT.input.contigs))
        | splitCsv(sep: '\t', header: true)
        | map { it -> [it.SAMPLE, file(it.CONTIGS)] }
        | set { contigs }

    reads = channel.empty()
    if (params.steps.binningONT.input.containsKey("reads")) {
        channel.from(file(params.steps.binningONT.input.reads))
            | splitCsv(sep: '\t', header: true)
            | map { it -> [it.SAMPLE, file(it.READS)] }
            | set { reads }
    }

    graph = channel.empty()
    if (params.steps.binningONT.input.containsKey("graph")) {
        channel.from(file(params.steps.binningONT.input.graph))
            | splitCsv(sep: '\t', header: true)
            | map { it -> [it.SAMPLE, file(it.GRAPH)] }
            | set { graph }
    }

    headerMapping = channel.empty()
    if (params.steps.binningONT.input.containsKey("headerMapping")) {
        channel.from(file(params.steps.binningONT.input.headerMapping))
            | splitCsv(sep: '\t', header: true)
            | map { it -> [it.SAMPLE, file(it.MAPPING)] }
            | set { headerMapping }
    }

    assemblyInfo = channel.empty()
    if (params.steps.binningONT.input.containsKey("assemblyInfo")) {
        channel.from(file(params.steps.binningONT.input.assemblyInfo))
            | splitCsv(sep: '\t', header: true)
            | map { it -> [it.SAMPLE, file(it.INFO)] }
            | set { assemblyInfo }
    }

    quality = channel.empty()
    if (params.steps.binningONT.input.containsKey("quality")) {
        channel.from(file(params.steps.binningONT.input.quality))
            | splitCsv(sep: '\t', header: true)
            | map { it -> [it.SAMPLE, it.QUALITY] }
            | set { quality }
    }

    wSaveSettingsList(reads | map { it -> it[SAMPLE_IDX] })

    _wBinning(contigs, reads, graph, headerMapping, assemblyInfo, quality)
}



/*
* The the format for every channel is described in corresponding file entrypoint.
*
*/
workflow wLongReadBinningList {
    take:
    contigs
    inputReads
    inputGraph
    headerMapping
    assemblyInfo
    medianQuality

    main:
    _wBinning(contigs, inputReads, inputGraph, headerMapping, assemblyInfo, medianQuality)

    emit:
    binsStats = _wBinning.out.binsStats
    bins = _wBinning.out.bins
    mapping = _wBinning.out.mapping
    notBinnedContigs = _wBinning.out.notBinnedContigs
    unmappedReads = _wBinning.out.unmappedReads
    contigCoverage = _wBinning.out.contigCoverage
}


workflow _wRunBinningTools {
    take:
    inputGraph
    contigs
    mappedReads
    headerMapping
    assemblyInfo
    medianQuality

    main:
    SAMPLE_IDX = 0

    inputGraph
        | join(contigs, by: SAMPLE_IDX)
        | join(mappedReads, by: SAMPLE_IDX)
        | join(headerMapping, by: SAMPLE_IDX)
        | join(assemblyInfo, by: SAMPLE_IDX)
        | pMetaCoAG

    pMetabat(
        channel.value(params?.steps?.containsKey("binningONT") && params?.steps?.binningONT.containsKey("metabat")),
        channel.value(
            [
                params.modules.binningONT,
                "metabatONT",
                params.steps?.binningONT?.metabat?.additionalParams,
            ]
        ),
        contigs | join(mappedReads, by: SAMPLE_IDX) | join(medianQuality, by: SAMPLE_IDX),
    )

    pSemiBin2(
        channel.value(params?.steps?.containsKey("binningONT") && params?.steps?.binningONT.containsKey("semibin2")),
        channel.value(
            [
                params.modules.binningONT,
                "semibin2ONT",
                params.steps?.binningONT?.semibin2?.additionalParams,
            ]
        ),
        contigs | join(mappedReads, by: SAMPLE_IDX),
    )

    pMetabat.out.bins
        | mix(pSemiBin2.out.bins)
        | mix(pMetaCoAG.out.bins)
        | set { bins }


    pMetabat.out.notBinned
        | mix(pSemiBin2.out.notBinned)
        | mix(pMetaCoAG.out.notBinned)
        | set { notBinned }

    // Compute bin statistcs (e.g. N50, average coverage depth, etc. ...)
    pMetabat.out.binContigMapping
        | join(mappedReads, by: SAMPLE_IDX)
        | combine(channel.from("metabatONT"))
        | join(pMetabat.out.bins, by: SAMPLE_IDX)
        | set { metabatBinStatisticsInput }
    pSemiBin2.out.binContigMapping
        | join(mappedReads, by: SAMPLE_IDX)
        | combine(channel.from("semibin2ONT"))
        | join(pSemiBin2.out.bins, by: SAMPLE_IDX)
        | set { semibin2BinStatisticsInput }
    pMetaCoAG.out.binContigMapping
        | join(mappedReads, by: SAMPLE_IDX)
        | combine(channel.from("metacoag"))
        | join(pMetaCoAG.out.bins, by: SAMPLE_IDX)
        | set { metacoagBinStatisticsInput }

    metabatBinStatisticsInput
        | mix(metacoagBinStatisticsInput)
        | mix(semibin2BinStatisticsInput)
        | join(medianQuality, by: SAMPLE_IDX)
        | set { binStatsInput }

    emit:
    bins = bins
    notBinned = notBinned
    binStatsInput = binStatsInput
}
/*
* The format of every channel is described in the corresponding File entrypoint.
*
*/
workflow _wBinning {
    take:
    contigs
    inputReads
    inputGraph
    headerMapping
    assemblyInfo
    medianQuality

    main:
    // Map reads against assembly and retrieve mapping quality
    SAMPLE_IDX = 0

    pMinimap2(
        channel.value(params?.steps?.containsKey("binningONT")),
        channel.value(
            [
                params.modules.binningONT,
                "contigMapping",
                params.steps?.binningONT?.minimap?.additionalParams?.minimap,
                params.steps?.binningONT?.minimap?.additionalParams?.samtoolsView,
                params.steps.containsKey("fragmentRecruitment"),
            ]
        ),
        contigs | join(inputReads, by: SAMPLE_IDX),
    )

    pMinimap2.out.mappedReads | set { mappedReads }
    wGetMappingQuality(params.modules.binningONT, mappedReads)

    pCovermContigsCoverage(
        channel.value(params?.steps?.binningONT.find { it.key == "contigsCoverage" }?.value),
        channel.value(
            [
                params.modules.binningONT,
                "contigCoverage",
                params?.steps?.binningONT?.contigsCoverage?.additionalParams,
            ]
        ),
        mappedReads | join(medianQuality, by: SAMPLE_IDX),
    )

    _wRunBinningTools(inputGraph, contigs, mappedReads, headerMapping, assemblyInfo, medianQuality)
    _wRunBinningTools.out.bins | set { bins }
    _wRunBinningTools.out.notBinned | set { notBinned }
    _wRunBinningTools.out.binStatsInput | set { binStatsInput }

    emptyFile = file(params.tempdir + "/empty")
    emptyFile.text = ""

    ALIGNMENT_INDEX = 2
    pCovermGenomeCoverage(
        channel.value(params?.steps?.binningONT.find { it.key == "genomeCoverage" }?.value),
        channel.value(""),
        channel.value(
            [
                params.modules.binningONT,
                "genomeCoverage",
                params?.steps?.binningONT?.genomeCoverage?.additionalParams,
            ]
        ),
        mappedReads | join(bins, by: SAMPLE_IDX) | map { sample ->
            sample.addAll(ALIGNMENT_INDEX, emptyFile)
            sample
        } | join(medianQuality, by: SAMPLE_IDX),
    )

    pCovermGenomeCoverage.out.logs | pDumpLogs

    // Flatten metabat outputs per sample and create a map with the 
    // following entries [BIN_ID:bin.name, SAMPLE:sample, PATH:bin]
    bins
        | map { it -> Utils.flattenTuple(it) }
        | flatMap { it -> createMap(it) }
        | set { binMap }

    pGetBinStatistics(channel.value(params.modules.binningONT), binStatsInput)

    // Add bin statistics 
    pGetBinStatistics.out.binsStats
        | map { _sample, stats -> file(stats) }
        | splitCsv(sep: '\t', header: true)
        | set { binsStats }
    mapJoin(binsStats, binMap, "BIN_ID", "BIN_ID") | set { binMap }

    emit:
    binsStats = binMap
    bins = bins
    mapping = mappedReads
    notBinnedContigs = notBinned
    unmappedReads = pMinimap2.out.unmappedReads
    contigCoverage = pCovermContigsCoverage.out.coverage
}
