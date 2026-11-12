include { pCollectFile; _wSplit } from './processes'
include { pDumpLogs } from '../utils/processes'

/**
*
* DIAMOND is used to search for big input queries in large protein databases.
* Multiple databases can be searched at the same time.
*
**/
process pDiamond {

      container "${params.diamond_image}"

      containerOptions Utils.constructParametersObject("diamond", params, params.apptainer) + (params.apptainer ? "" : " --entrypoint='' " + Utils.getDockerNetwork())

      memory { Utils.getMemoryResources(params.resources.large, "${sample}", task.attempt, params.resources) }

      cpus { Utils.getCPUsResources(params.resources.large, "${sample}", task.attempt, params.resources) }

      tag "Sample: $sample, Database: $dbType"

      secret { "${S3_DB_ACCESS}" != "" ? ["S3_${dbType}_ACCESS", "S3_${dbType}_SECRET"] : [] }

   input:
      val(binType)
      tuple val(sample), file(fasta), file(contig2GeneMapping), val(start), val(stop), val(numberOfChunks), val(dbType), val(parameters), val(columns), val(EXTRACTED_DB), val(DOWNLOAD_LINK), val(MD5SUM), val(S5CMD_PARAMS)

   output:
      tuple val("${dbType}"), val("${sample}"), val("${binType}"), val("${start}"), val("${stop}"), val("${numberOfChunks}"),\
        path("${sample}_${binType}.${start}.${stop}.${dbType}.blast.tsv"), optional: true, emit: blast
      tuple val("${sample}_${binType}"), val("${output}"), val(params.LOG_LEVELS.INFO), file(".command.sh"), \
        file(".command.out"), file(".command.err"), file(".command.log"), emit: logs

   script:
   output = Output.getOutput("${sample}", params.runid, "diamond/${dbType}", params.modules.annotation, "")
   S3_DB_ACCESS = params.steps?.annotation?.diamond?."${dbType}"?.database?.download?.s5cmd && S5CMD_PARAMS.indexOf("--no-sign-request") == -1 ? "\$S3_${dbType}_ACCESS" : ""
   S3_DB_SECRET = params.steps?.annotation?.diamond?."${dbType}"?.database?.download?.s5cmd && S5CMD_PARAMS.indexOf("--no-sign-request") == -1 ? "\$S3_${dbType}_SECRET" : ""
   ADDITIONAL_COLUMNS = "${columns}".trim() == "" ? "" : ",${columns}"
   template("diamond.sh")
}

/**
*
* It searches the predicted proteins against all the configured databases in the params.steps.annotation.diamond section.
* The proteins are split into chunks and searched in parallel. The results of these chunks are then collected and split into one file per bin.
*
* sourceChannel: The type of the input (e.g. binned, unbinned).
* contig2GeneMapping: Mapping of locus tags to contigs and bins per sample.
* proteins: Protein fasta files per sample and bin.
* fastaCounter: The expected number of bins or fasta files per sample.
*
**/
workflow _wAnalyseDiamond {
  take:
    sourceChannel
    contig2GeneMapping
    proteins
    fastaCounter
  main:
    SAMPLE_IDX = 0
    PATH_IDX = 2

    DIAMOND_CHUNK_SIZE_DEFAULT = 7000

    // Collect all databases
    selectedDBs = (params?.steps?.annotation?.diamond ?: [:])
        .findAll { it.key != "chunkSize" }
        .collect {
            [
                it.key,
                it.value?.additionalParams?.search ?: "",
                it.value?.additionalParams?.additionalColumns ?: "",
                it.value?.database?.extractedDBPath ?: "",
                it.value?.database?.download?.source ?: "",
                it.value?.database?.download?.md5sum ?: "",
                it.value?.database?.download?.s5cmd?.params ?: "",
            ]
        }

    // The following groupTuple operators need a counter for the expected number of bins or fasta files of a sample.
    // This allows groupTuple to not block until all samples are processed.

    // Run all amino acid outputs against all databases
    // Collect by sample name to bundle searches and avoid calls with small input files
    proteins
        | map { [it[SAMPLE_IDX], it[PATH_IDX]] }
        | combine(fastaCounter, by: SAMPLE_IDX)
        | map { sample, path, size -> tuple(groupKey(sample, size), path) }
        | groupTuple()
        | combine(contig2GeneMapping, by: SAMPLE_IDX)
        | set { groupedProteins }

    // Split fasta files and run DIAMOND on each part in parallel
    _wSplit(groupedProteins, Channel.from(params.steps?.annotation?.diamond?.chunkSize ?: DIAMOND_CHUNK_SIZE_DEFAULT))
        | combine(Channel.from(selectedDBs))
        | set { combinedDiamond }

    pDiamond(sourceChannel, combinedDiamond)

    FIRST_ELEM_IDX = 0
    FIRST_ELEM_DB_IDX = 0
    FIRST_ELEM_SAMPLE_IDX = 1
    FIRST_ELEM_TYPE_IDX = 2

    ELEM_PATH_IDX = 3

    pDiamond.out.blast
        // An artificial group key is created which consists of the db, sample name and the type of the data (e.g. binned)
        | map { db, sample, type, start, stop, chunks, out -> [db + "_-_" + sample + "_-_" + type, db, sample, type, start, stop, chunks, out] }
        | map { key, db, sample, type, start, stop, chunks, out -> tuple(groupKey(key, chunks.toInteger()), [db, sample, type, out]) }
        | groupTuple()
        | map { key, dataset ->
            [
                dataset[FIRST_ELEM_IDX][FIRST_ELEM_DB_IDX],
                dataset[FIRST_ELEM_IDX][FIRST_ELEM_SAMPLE_IDX],
                dataset[FIRST_ELEM_IDX][FIRST_ELEM_TYPE_IDX],
                dataset.stream().map { elem -> elem[ELEM_PATH_IDX] }.collect(),
            ]
        }
        | combine(Channel.value("diamond"))
        | combine(Channel.value("blast"))
        | pCollectFile

    pCollectFile.out.blast
        | flatMap { sample, type, dbType, blastFiles -> blastFiles.stream().map { fi -> [sample, type, dbType, fi] }.collect() }
        | set { collectedDiamondResults }

    pDiamond.out.logs | pDumpLogs

  emit:
    blast = collectedDiamondResults
}
