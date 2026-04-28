include {
    pGetBinStatistics ;
    pCovermContigsCoverage ;
    pCovermGenomeCoverage ;
    pSemiBin2 ;
    pMetabat ;
    _wRunMappers
} from './processes'
include { pProdigal ; pHmmSearch } from '../annotation/module'
include { wSaveSettingsList } from '../config/module'
include { createMap ; mapJoin } from '../utils/methods'

include { pDumpLogs } from '../utils/processes'


process pMetabinner {

    container "${params.metabinner_image}"

    tag "Sample: ${sample}"

    label 'highmemLarge'

    publishDir params.output, mode: "${params.publishDirMode}", saveAs: { filename ->
        Output.getOutput("${sample}", params.runid, "metabinner", params.modules.binning, filename)
    }

    when params.steps.containsKey("binning") && params.steps.binning.containsKey("metabinner")

    containerOptions params.apptainer ? "" : ' --user 0:0 '

    input:
    tuple val(sample), path(contigs), path(bam)

    output:
    tuple val("${sample}"), path("${sample}_bin.*.fa", arity: '1..*'), optional: true, emit: bins
    tuple val("${sample}"), file("${sample}_notBinned.fa"), optional: true, emit: notBinned
    tuple val("${sample}"), file("${sample}_bin_contig_mapping.tsv"), optional: true, emit: binContigMapping
    tuple file(".command.sh"), file(".command.out"), file(".command.err"), file(".command.log")

    shell:
    template('metabinner.sh')
}


/**
 * MAGScoT - Run MAGScoT binning refinement
 * @param sample: Sample name
 * @param contigMaps: Individual contig to bin files from binning algorithms
 * @param allHits: GTDBtk single copy marker identification of all contigs
 * @param contigs: Path to contigs
 * @return: Files containing scores, bin-contig mappings, new bins, and not binned contigs
 *
 * This process runs the MAGScoT R script on the input contig maps
 * to compare and refine the binning of multiple other binners.
 * It outputs several files containing scores, bin-contig mappings, bins, and not binned contigs.
 **/
process pMAGScoT {

    container "${params.magscot_image}"

    tag "Sample: ${sample}"

    label 'small'

    publishDir params.output, mode: "${params.publishDirMode}", saveAs: { filename ->
        Output.getOutput("${sample}", params.runid, "magscot", params.modules.binning, filename)
    }

    // Override default MAGScoT container entrypoint, so that the "RScript magscot" call that is normally run
    // does not clash with the Nextflow process call "/bin/bash"
    containerOptions params.apptainer ? "" : ' --entrypoint "" '

    when params.steps.containsKey("binning") && params.steps.binning.containsKey("magscot")

    input:
    tuple val(sample), file(contigMaps), file(allHits), path(contigs)

    output:
    tuple val("${sample}"), file("${sample}_MagScoT.*"), optional: true, emit: scores
    tuple val("${sample}"), file("${sample}_bin_contig_mapping.tsv"), optional: true, emit: binContigMapping
    tuple val("${sample}"), path("${sample}_bin.*.fa", arity: '1..*'), optional: true, emit: bins
    tuple val("${sample}"), file("${sample}_notBinned.fa"), optional: true, emit: notBinned
    tuple file(".command.sh"), file(".command.out"), file(".command.err"), file(".command.log")

    shell:
    '''
    # Once in a blue moon Nextflow leaves the header in the file
    # Failsafe to remove header from contigMaps file if it exists
    sed -i '/^BIN_ID\tCONTIG\tBINNER$/d' !{contigMaps}
    Rscript /opt/MAGScoT.R !{params.steps?.binning?.magscot?.additionalParams} -i !{contigMaps} --hmm !{allHits} -o !{sample}_MagScoT

    # Create a new binning file according to the naming convention
    echo "Converting MAGScoT binning according to the naming convention"
    echo -e 'BIN_ID\tCONTIG\tBINNER' > !{sample}_bin_contig_mapping.tsv

    # Use head to get the first line, awk to print the first field (the filename),
    # and sed to remove everything up to and including the last dot (.) in the filename, leaving only the file extension.
    # Export the variable to use it in the xargs commands subshell further down
    export EXT=$(head -n 1 !{contigMaps} | awk '{print $1}' | sed 's/.*\\.//')

    # Remove the header and pipe the remaining lines to xargs to run the script line by line
    sed 1d !{sample}_MagScoT.refined.contig_to_bin.out | xargs -n 2 sh -c '
        # Get the first column and separate the number, remove the leading zeros
        binID=$(echo $0 | cut -d"_" -f4 | sed 's/^0*//')
        CONTIG=$1
        # Create a file with the contigs for each bin for reconstruction with seqkit
        echo $CONTIG >> !{sample}_bin.$binID.lst
        echo "!{sample}_bin.$binID.$EXT\t$CONTIG\tMAGScot" >> !{sample}_bin_contig_mapping.tsv
    '
    echo "Done converting"
    # Reconstructing the bins from the converted MAGScoT output list
    echo "Reconstructing bins"
    for bin in !{sample}_bin.*.lst; do
        seqkit grep -f $bin !{contigs} -o ${bin%.lst}.fa.tmp
    done
    echo "Done reconstructing bins"

    # Rename the bins to the naming convention
    for bin in $(find * -name "*bin*.fa.tmp"); do
    	BIN_NAME="$(basename ${bin})"
    	echo "Renaming bin: ${BIN_NAME}"

    	# Get id of the bin (e.g get 2 of the bin SAMPLEID_bin.2.fa.tmp)
    	ID=$(echo ${BIN_NAME} | rev | cut -d '.' -f 3 | rev)
    	echo "Bin id: ${ID}"

    	# Append bin id to every header
    	seqkit replace  -p '(.*)' -r "\\${1} MAG=${ID}" $bin > ${BIN_NAME%.tmp}
    	# Remove temporary file
    	rm $bin
    done

    # Creating un-binned contigs file
    seqkit grep -f <(cat !{sample}_bin.*.lst) -v !{contigs} -o !{sample}_notBinned.fa.tmp
    seqkit replace -p '(.*)' -r "\\${1} MAG=NotBinned" !{sample}_notBinned.fa.tmp > !{sample}_notBinned.fa
    '''
}


/*
 * Takes three tab separated file of files containing contigs, paired and optional single reads 
 * as input and produces binning results.
 * Input files for reads must have two columns seperated by tabs:
 * SAMPLE and READS
 * Input files for contigs must have two columns seperated by tabs:
 * SAMPLE and CONTIGS
 *
 */
workflow wShortReadBinningFile {
    SAMPLE_IDX = 0
    SAMPLE_PAIRED_IDX = 1
    UNPAIRED_IDX = 2

    channel.from(file(params.steps.binning.input.contigs))
        | splitCsv(sep: '\t', header: true)
        | map { it -> [it.SAMPLE, file(it.CONTIGS)] }
        | set { contigs }

    readsPaired = channel.empty()
    if (params.steps.binning.input.containsKey("paired")) {
        channel.from(file(params.steps.binning.input.paired))
            | splitCsv(sep: '\t', header: true)
            | map { it -> [it.SAMPLE, file(it.READS)] }
            | set { readsPaired }
    }

    readsSingle = channel.empty()
    if (params.steps.binning.input.containsKey("single")) {
        channel.from(file(params.steps.binning.input.single))
            | splitCsv(sep: '\t', header: true)
            | map { it -> [it.SAMPLE, file(it.READS)] }
            | set { readsSingle }
    }

    readsPaired
        | join(readsSingle, by: SAMPLE_IDX, remainder: true)
        | map { sample ->
            sample[UNPAIRED_IDX] == null
                ? [sample[SAMPLE_IDX], sample[SAMPLE_PAIRED_IDX], file("NOT_SET")]
                : sample
        }
        | set { reads }

    wSaveSettingsList(reads | map { it -> it[SAMPLE_IDX] })

    _wBinning(contigs, reads)
}


/*
* This workflow accepts channels containing contigs and reads with the following format:
*  
* Reads: SAMPLE and READS
* Contigs: SAMPLE and CONTIGS
*
*/
workflow wShortReadBinningList {
    take:
    contigs
    inputReads

    main:
    _wBinning(contigs, inputReads)

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
    contigs
    mappedReads

    main:
    DO_NOT_ESTIMATE_IDENTITY = "-1"
    SAMPLE_IDX = 0

    contigs | join(mappedReads, by: SAMPLE_IDX) | set { binningInput }
    pMetabinner(binningInput)

    pSemiBin2(
        channel.value(params?.steps?.containsKey("binning") && params?.steps?.binning.containsKey("semibin2")),
        channel.value(
            [
                params.modules.binning,
                "semibin2",
                params.steps?.binning?.semibin2?.additionalParams,
            ]
        ),
        binningInput,
    )

    pMetabat(
        channel.value(params?.steps?.containsKey("binning") && params?.steps?.binning.containsKey("metabat")),
        channel.value(
            [
                params.modules.binning,
                "metabat",
                params.steps?.binning?.metabat?.additionalParams,
            ]
        ),
        binningInput | combine(channel.value(DO_NOT_ESTIMATE_IDENTITY)),
    )

    pMetabinner.out.bins
        | mix(pMetabat.out.bins)
        | mix(pSemiBin2.out.bins)
        | set { bins }
    pMetabinner.out.notBinned
        | mix(pSemiBin2.out.notBinned)
        | mix(pMetabat.out.notBinned)
        | set { notBinned }

    pMetabinner.out.binContigMapping
        | join(mappedReads, by: SAMPLE_IDX)
        | combine(channel.from("metabinner"))
        | join(pMetabinner.out.bins, by: SAMPLE_IDX)
        | set { metabinnerBinStatisticsInput }
    pMetabat.out.binContigMapping
        | join(mappedReads, by: SAMPLE_IDX)
        | combine(channel.from("metabat"))
        | join(pMetabat.out.bins, by: SAMPLE_IDX)
        | set { metabatBinStatisticsInput }
    pSemiBin2.out.binContigMapping
        | join(mappedReads, by: SAMPLE_IDX)
        | combine(channel.from("semibin2"))
        | join(pSemiBin2.out.bins, by: SAMPLE_IDX)
        | set { semibin2BinStatisticsInput }

    metabatBinStatisticsInput
        | mix(semibin2BinStatisticsInput)
        | mix(metabinnerBinStatisticsInput)
        | combine(channel.value(DO_NOT_ESTIMATE_IDENTITY))
        | set { binStatsInput }

    emit:
    bins = bins
    notBinned = notBinned
    binStatsInput = binStatsInput
}


/*
*
* This workflow takes an input_reads channel as input with the following format [SAMPLE, READS PAIRED, READS UNPAIRED]
* and a contigs channel with the format [SAMPLE, CONTIGS]
*
*/
workflow _wBinning {
    take:
    contigs
    inputReads

    main:
    // Map reads against assembly and retrieve mapping quality
    SAMPLE_IDX = 0
    DO_NOT_ESTIMATE_IDENTITY = "-1"

    mappedReads = channel.empty()
    unmappedReads = channel.empty()
    if (params.steps.containsKey("binning")) {
        def activeMapper = params.steps?.binning?.keySet()?.find { tool -> ["bwa", "bowtie", "bwa2"].contains(tool) }
        def mapper = params.steps?.binning[activeMapper]

        mapperConfig = mapper
            ? channel.value(
                [
                    params.modules.binning,
                    "contigMapping",
                    mapper.additionalParams[activeMapper],
                    mapper.additionalParams.samtoolsView,
                    params.steps.containsKey("fragmentRecruitment"),
                ]
            )
            : channel.empty()

        contigs | combine(inputReads, by: SAMPLE_IDX) | set { mapperInput }

        _wRunMappers(params.modules.binning, channel.value(activeMapper), mapperConfig, mapperInput)

        _wRunMappers.out.mappedReads | set { mappedReads }

        _wRunMappers.out.unmappedReads | set { unmappedReads }
    }


    pCovermContigsCoverage(
        channel.value(params?.steps?.binning.find { it.key == "contigsCoverage" }?.value),
        channel.value(
            [
                params.modules.binning,
                "contigCoverage",
                params?.steps?.binning?.contigsCoverage?.additionalParams,
            ]
        ),
        mappedReads | combine(channel.value(DO_NOT_ESTIMATE_IDENTITY)),
    )

    // RUn the actual Binning tools
    _wRunBinningTools(contigs, mappedReads)
    _wRunBinningTools.out.bins | set { bins }
    _wRunBinningTools.out.notBinned | set { notBinned }
    _wRunBinningTools.out.binStatsInput | set { binStatsInput }

    // Re-evaluate binning with MAGScoT
    // ORF detection with Prodigal for MAGScoT
    CONTIG_MAPPING_IDX = 1

    magscot_input = channel.empty()
    if (params.steps.containsKey("binning") && params.steps.binning.containsKey("magscot")) {
        pProdigal(contigs)
        pHmmSearch(pProdigal.out.prodigal_faa)
        pMetabinner.out.binContigMapping
            | mix(pSemiBin2.out.binContigMapping)
            | mix(pMetabat.out.binContigMapping)
            | collectFile(keepHeader: false) { item -> ["${item[SAMPLE_IDX]}", item[CONTIG_MAPPING_IDX].text] }
            | map { f -> [file(f).name, f] }
            | join(pHmmSearch.out.allhits, by: SAMPLE_IDX)
            | join(contigs, by: SAMPLE_IDX)
            | set { magscot_input }
    }
    pMAGScoT(magscot_input)

    // Only use MAGScoT bins if the user has selected the refinement step
    if (params.steps.containsKey("binning") && params.steps.binning.containsKey("magscot")) {
        pMAGScoT.out.bins | set { bins }
        pMAGScoT.out.notBinned | set { notBinned }
    }

    emptyFile = file(params.tempdir + "/empty")

    ALIGNMENT_INDEX = 2
    pCovermGenomeCoverage(
        channel.value(params?.steps?.binning.find { it.key == "genomeCoverage" }?.value),
        channel.value(""),
        channel.value(
            [
                params.modules.binning,
                "genomeCoverage",
                params?.steps?.binning?.genomeCoverage?.additionalParams,
            ]
        ),
        mappedReads | join(bins, by: SAMPLE_IDX) | map { sample ->
            sample.addAll(ALIGNMENT_INDEX, emptyFile)
            sample
        } | combine(channel.value(DO_NOT_ESTIMATE_IDENTITY)),
    )

    pCovermGenomeCoverage.out.logs | pDumpLogs

    // Flatten binning outputs per sample and create a map with the 
    // following entries [BIN_ID:bin.name, SAMPLE:sample, PATH:bin]
    bins | map { it -> Utils.flattenTuple(it) } | flatMap { it -> createMap(it) } | set { binMap }

    // Compute bin statistcs (e.g. N50, average coverage depth, etc. ...)
    // Only use the MAGScoT bin statistics if the user has selected the refinement step
    if (params.steps.containsKey("binning") && params.steps.binning.containsKey("magscot")) {
        pMAGScoT.out.binContigMapping
            | join(mappedReads, by: SAMPLE_IDX)
            | combine(channel.from("magscot"))
            | join(pMAGScoT.out.bins, by: SAMPLE_IDX)
            | set { magscotBinStatisticsInput }

        magscotBinStatisticsInput | combine(channel.value(DO_NOT_ESTIMATE_IDENTITY)) | set { binStatsInput }
    }

    pGetBinStatistics(channel.value(params.modules.binning), binStatsInput)

    // Add bin statistics 
    pGetBinStatistics.out.binsStats
        | map { it -> file(it[1]) }
        | splitCsv(sep: '\t', header: true)
        | set { binsStats }
    mapJoin(binsStats, binMap, "BIN_ID", "BIN_ID") | set { binMap }

    emit:
    binsStats = binMap
    bins = bins
    mapping = mappedReads
    notBinnedContigs = notBinned
    unmappedReads = unmappedReads
    contigCoverage = pCovermContigsCoverage.out.coverage
}
