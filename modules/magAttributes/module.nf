include { wSaveSettingsList } from '../config/module'

include { _wCollectFile as _wCollectCheckm; \
        _wCollectFile as _wCollectCheckm2; \
        _wCollectFile as _wCollectGtdbtkArchea
	_wCollectFile as _wCollectGtdbtkBac
	_wCollectFile as _wCollectGtdbtkUnclassified
	_wCollectFile as _wCollectGtdbtkCombined
	_wCollectFile as _wCollectGtdbtkMissing
	_wCollectFile as _wCollectGtdbtkSummary; } from './processes'
include { pDumpLogs } from '../utils/processes'
include { pProkka; _wCreateProkkaInput } from '../annotation/module'

def getOutput(SAMPLE, RUNID, TOOL, filename){
    return SAMPLE + '/' + RUNID + '/' + params.modules.magAttributes.name + '/' + 
         params.modules.magAttributes.version.major + "." +  
         params.modules.magAttributes.version.minor + "." +
         params.modules.magAttributes.version.patch +
         '/' + TOOL + '/' + filename
}


process pCmseq {

    container "${params.cmseq_image}"

    label 'tiny'

    when params.steps.magAttributes.containsKey("cmseq")

    input:
    tuple val(sample), val(bin), file(gff), file(bam), file(bai)

    output:
    file("${sample}_${bin}.txt")

    shell:
    '''
    zcat -f !{gff} > input.gff
    polymut.py --mincov 10 --gff_file input.gff !{bam} > !{sample}_!{bin}.txt
    '''
}


process pCheckM {

    container "${params.checkm_image}"

    tag "Sample: $sample"

    secret { "${S3_checkm_ACCESS}"!="" ? ["S3_checkm_ACCESS", "S3_checkm_SECRET"] : [] } 

    publishDir params.output, mode: "${params.publishDirMode}", saveAs: { filename -> getOutput("${sample}", params.runid, "checkm", filename) }, \
      pattern: "{}"

    when params.steps.containsKey("magAttributes") && params.steps.magAttributes.containsKey("checkm") \
	&& !params.steps.magAttributes.containsKey("checkm2")

    containerOptions Utils.getDockerMount(params.steps?.magAttributes?.checkm?.database, params, apptainer=params.apptainer) + (params.apptainer ? "" : Utils.getDockerNetwork())

    beforeScript Utils.getCreateDatabaseDirCommand("${params.polished.databases}")

    label 'highmemMedium'

    input:
    tuple val(sample), val(ending), path(bins), val(chunkId), val(type), val(size)

    output:
    tuple val("${sample}"), path("${sample}_checkm_*.tsv", type: "file"), val(size), emit: checkm
    tuple env(FILE_ID), val("${output}"), val(params.LOG_LEVELS.INFO), file(".command.sh"), \
	file(".command.out"), file(".command.err"), file(".command.log"), emit: logs

    shell:
    output = getOutput("${sample}", params.runid, "checkm", "")
    S5CMD_PARAMS=params?.steps?.magAttributes?.checkm?.database?.download?.s5cmd?.params ?: "" 
    DOWNLOAD_LINK=params?.steps?.magAttributes?.checkm?.database?.download?.source ?: ""
    MD5SUM=params.steps?.magAttributes?.checkm?.database?.download?.md5sum ?: ""
    EXTRACTED_DB=params.steps?.magAttributes?.checkm?.database?.extractedDBPath ?: ""
    S3_checkm_ACCESS=params?.steps?.magAttributes?.checkm?.database?.download?.s5cmd && S5CMD_PARAMS.indexOf("--no-sign-request") == -1 ? "\$S3_checkm_ACCESS" : ""
    S3_checkm_SECRET=params?.steps?.magAttributes?.checkm?.database?.download?.s5cmd && S5CMD_PARAMS.indexOf("--no-sign-request") == -1 ? "\$S3_checkm_SECRET" : ""
    FILE_TYPE = type == "" ? "" : type + "_"
    template 'checkm.sh'
}


process pCheckM2 {

    container "${params.checkm2_image}"

    tag "Sample: $sample"

    secret { "${S3_checkm2_ACCESS}"!="" ? ["S3_checkm2_ACCESS", "S3_checkm2_SECRET"] : [] } 

    publishDir params.output, mode: "${params.publishDirMode}", saveAs: { filename -> getOutput("${sample}", params.runid, "checkm2", filename) }, \
      pattern: "{}"

    when params.steps.containsKey("magAttributes") && params.steps.magAttributes.containsKey("checkm2")

    containerOptions  Utils.getDockerMount(params.steps?.magAttributes?.checkm2?.database, params, apptainer=params.apptainer) + (params.apptainer ? "" : Utils.getDockerNetwork())

    beforeScript Utils.getCreateDatabaseDirCommand("${params.polished.databases}")

    label 'medium'

    input:
    tuple val(sample), val(ending), path(bins), val(chunkId), val(type), val(size)

    output:
    tuple val("${sample}"), path("${sample}_checkm2_*.tsv", type: "file"), val(size), emit: checkm
    tuple env(FILE_ID), val("${output}"), val(params.LOG_LEVELS.INFO), file(".command.sh"), \
	file(".command.out"), file(".command.err"), file(".command.log"), emit: logs

    shell:
    output = getOutput("${sample}", params.runid, "checkm2", "")
    S5CMD_PARAMS=params?.steps?.magAttributes?.checkm2?.database?.download?.s5cmd?.params ?: "" 
    DOWNLOAD_LINK=params?.steps?.magAttributes?.checkm2?.database?.download?.source ?: ""
    MD5SUM=params.steps?.magAttributes?.checkm2?.database?.download?.md5sum ?: ""
    EXTRACTED_DB=params.steps?.magAttributes?.checkm2?.database?.extractedDBPath ?: ""
    S3_checkm2_ACCESS=params?.steps?.magAttributes?.checkm2?.database?.download?.s5cmd && S5CMD_PARAMS.indexOf("--no-sign-request") == -1 ? "\$S3_checkm2_ACCESS" : ""
    S3_checkm2_SECRET=params?.steps?.magAttributes?.checkm2?.database?.download?.s5cmd && S5CMD_PARAMS.indexOf("--no-sign-request") == -1 ? "\$S3_checkm2_SECRET" : ""
    FILE_TYPE = type == "" ? "" : type + "_"
    template 'checkm2.sh'
}

process pGtdbtk {

    container "${params.gtdbtk_image}"

    label 'highmemMedium'

    tag "Sample: $sample"

    secret { "${S3_gtdb_ACCESS}"!="" ? ["S3_gtdb_ACCESS", "S3_gtdb_SECRET"] : [] } 

    publishDir params.output, mode: "${params.publishDirMode}", saveAs: { filename -> getOutput("${sample}",params.runid ,"gtdb", filename) }, \
      pattern: "{**.tree}"

    when params.steps.containsKey("magAttributes") && params.steps.magAttributes.containsKey("gtdb")

    containerOptions Utils.getDockerMount(params?.steps?.magAttributes?.gtdb?.database, params, apptainer=params.apptainer) + (params.apptainer ? "" : Utils.getDockerNetwork())

    beforeScript Utils.getCreateDatabaseDirCommand("${params.polished.databases}")

    input:
    tuple val(sample), val(ending), path(bins), val(chunkId), val(type), val(size)

    output:
    tuple val("${sample}"), path("chunk_*_${sample}_gtdbtk.bac120.summary.tsv"), val(size), optional: true, emit: bacteria
    tuple val("${sample}"), path("chunk_*_${sample}_gtdbtk.ar122.summary.tsv"), val(size), optional: true, emit: archea
    tuple val("${sample}"), path("chunk_*_${sample}_gtdbtk_summary_raw_combined.tsv"), val(size), optional: true, emit: summaryRawCombined
    tuple val("${sample}"), path("chunk_*_${sample}_gtdbtk_unclassified.tsv"), val(size), optional: true, emit: unclassified
    tuple val("${sample}"), path("*.tree"), val("${sample}"), val(size), optional: true, emit: tree
    tuple val("${sample}"), path("chunk_*_${sample}_gtdbtk_combined.tsv"), val(size), optional: true, emit: combined
    tuple val("${sample}"), path("chunk_*_${sample}_missing_bins.tsv"), val(size), optional: true, emit: missing
    tuple env(FILE_ID), val("${output}"), val(params.LOG_LEVELS.INFO), file(".command.sh"), \
	file(".command.out"), file(".command.err"), file(".command.log"), emit: logs

    shell:
    output = getOutput("${sample}", params.runid, "gtdb", "")
    S5CMD_PARAMS=params?.steps?.magAttributes?.gtdb?.database?.download?.s5cmd?.params ?: ""
    DOWNLOAD_LINK=params?.steps?.magAttributes?.gtdb?.database?.download?.source ?: ""
    MD5SUM=params.steps?.magAttributes?.gtdb?.database?.download?.md5sum ?: ""
    EXTRACTED_DB=params.steps?.magAttributes?.gtdb?.database?.extractedDBPath ?: ""
    GTDB_PARAMS=params.steps.magAttributes.gtdb.additionalParams ?: ""
    S3_gtdb_ACCESS=params?.steps?.magAttributes?.gtdb?.database?.download?.s5cmd && S5CMD_PARAMS.indexOf("--no-sign-request") == -1 ? "\$S3_gtdb_ACCESS" : ""
    S3_gtdb_SECRET=params?.steps?.magAttributes?.gtdb?.database?.download?.s5cmd && S5CMD_PARAMS.indexOf("--no-sign-request") == -1 ? "\$S3_gtdb_SECRET" : ""
    FILE_TYPE = type == "" ? "" : type + "_"
    template 'gtdb.sh'
}


/*
* The CMSeq workflow should estimate the heterogenety of a MAG in a given sample.
* Input:
*     * genomes  - A table with the columns PATH and  DATASET
*         Example:
*          PATH    DATASET
*          /vol/spool/cmseq_20210607/GCF_000145295.fa      SAMPLE
* 
*     * alignments - A table with read alignment against genomes. 
*           Example:
*           BAM	PATH    BAI
*           SAMPLE /vol/spool/fragmentRecruitment_20210607/work/3d/9a45b85c15b22a6bb8ed4635391a40/ERR2019981.bam   /vol/spool/fragmentRecruitment_20210607/work/3d/9a45b85c15b22a6bb8ed4635391a40/ERR2019981.bam.bai
*
*/
workflow wCMSeqWorkflowFile {
   take:
      genomes
      alignments
   main:
      wMagAttributesFile(genomes)

      
      //prepare alignment inputs for CMSeq:
      alignments | splitCsv(sep: '\t', header: true)  \
       | map { sample -> [sample.BAM, file(sample.PATH), file(sample.BAI)] } | set {alignments}

      //prepare genome inputs for prokka:
      genomes | splitCsv(sep: '\t', header: true) \
       | map { sample -> [sample.DATASET, file(sample.PATH).getName(), file(sample.PATH)] } \
       | set { genomesForProkka }

      DATASET_IDX = 0 
      wSaveSettingsList(samplesContigs | map { it -> it[DATASET_IDX] })

      _wCreateProkkaInput(genomesForProkka, wMagAttributesFile.out.gtdb)

      pProkka(Channel.value("meta"), _wCreateProkkaInput.out.prokkaInput)

      SAMPLE_IDX = 0
      SAMPLE_BIN_IDX = 1
      GFF_IDX = 2
      ALIGNMENT_IDX = 3
      ALIGNMENT_INDEX_IDX = 4

      pProkka.out.gff \
         | join(alignments, by: SAMPLE_IDX, remainder: true)   \
         | map { it -> [it[SAMPLE_IDX], it[SAMPLE_BIN_IDX], it[GFF_IDX], it[ALIGNMENT_IDX], it[ALIGNMENT_INDEX_IDX]] } |  pCmseq

}


/**
* This workflow provides the same functionality as wMagAttributesList.
* Input:
* It accepts as input all a tsv table of the following format:
* 
* DATASET	PATH	BIN_ID
* SAMPLE_NAME	/path/to/bin.fa	sample_name_bin
*
*/
workflow wMagAttributesFile {
   take:
     mags
   main:
     DATASET_IDX = 0
     mags | splitCsv(sep: '\t', header: true) \
       | map { sample -> [sample.DATASET, file(sample.PATH)] } \
       | groupTuple(by: DATASET_IDX) | set { mags }
  
     _wMagAttributes(Channel.value(""), mags)
     wSaveSettingsList(mags | map { it -> it[DATASET_IDX] })
   emit:
     checkm = _wMagAttributes.out.checkm
     gtdb = _wMagAttributes.out.gtdb
}


/**
*
* This module infers any general bin specific properties such as contamination, completeness, lineage.
*
* Input:
* It accepts a list where each entry has the following format: [SAMPLE_NAME, ["/path/to/bin.fa", "/path/to/bin.fa"]] 
*
* Output:
*  * Map of checkm values 
*    Example: 
*    [PATH:/vol/spool/meta/test/bins/small/bin.1.fa, SAMPLE:test1,
*     BIN_ID:bin.1, Marker lineage:k__Bacteria (UID203), # genomes:5449, # markers:103, # marker sets:57,
*     0:97, 1:6, 2:0, 3:0, 4:0, 5+:0, COMPLETENESS:7.72, CONTAMINATION:0.00, HETEROGENEITY:0.00] 
*
*/
workflow wMagAttributesList {
   take: 
     type
     mags
   main:
     _wMagAttributes(type, mags)
   emit:
     checkm = _wMagAttributes.out.checkm
     gtdb = _wMagAttributes.out.gtdb
     gtdbMissing = _wMagAttributes.out.gtdbMissing
     checkmFiles = _wMagAttributes.out.checkmFiles
     gtdbCombinedSummaryFiles = _wMagAttributes.out.gtdbCombinedSummaryFiles
}

/*
*
* Method takes a list of the form [SAMPLE, [BIN1 path, BIN2 path]] as input
* and produces a flattend list which is grouped by dataset and sample.
* The output has the form [SAMPLE, file ending (e.g. .fa), [BIN 1 path, BIN 2 path], chunk id]
*
*/
def groupBins(binning, buffer){
  def chunkList = [];
  def SAMPLE_IDX = 0;
  def BIN_PATHS_IDX = 1;
  def BIN_TYPE = 2;
  binning[BIN_PATHS_IDX].collate(buffer).eachWithIndex {
       it, indexSample -> it.groupBy{ 
            bin -> file(bin).name.substring(file(bin).name.lastIndexOf(".")) 
       }.eachWithIndex {
            ending, group, indexFileEnding -> \
	chunkList.add([binning[SAMPLE_IDX], ending, group, \
	indexSample.toString() + indexFileEnding.toString(), \
	binning[BIN_TYPE]]);
       }
  }
  
  for (def List chunk: chunkList) {
      chunk.add(chunkList.size())
  } 

  return chunkList;
}

workflow _wMagAttributes {
   take: 
     type
     bins
   main:
     GTDB_DEFAULT_BUFFER = 500
     CHECKM_DEFAULT_BUFFER = 30
     CHECKM2_DEFAULT_BUFFER = 20000
     BIN_FILES_INPUT_IDX = 1

     DATASET_IDX = 0
     FILE_ENDING_IDX = 1
     BIN_FILES_IDX = 2
     BIN_FILES_OUTPUT_IDX = 0

     bins | combine(type) | set { binsWithType } 

     // get file ending of bin files (.fa, .fasta, ...) and group by file ending and dataset
     binsWithType | flatMap({n -> groupBins(n, params?.steps?.magAttributes?.checkm?.buffer ?: CHECKM_DEFAULT_BUFFER)}) \
       | pCheckM | set { checkm }
     binsWithType | flatMap({n -> groupBins(n, params?.steps?.magAttributes?.checkm2?.buffer ?: CHECKM2_DEFAULT_BUFFER)}) \
       | pCheckM2 | set { checkm2 }
     binsWithType | flatMap({n -> groupBins(n, params?.steps?.magAttributes?.gtdb?.buffer ?: GTDB_DEFAULT_BUFFER )}) \
       | pGtdbtk | set { gtdb }

     toolCheckm = Channel.value("checkm")
     toolCheckm2 = Channel.value("checkm2")
     _wCollectCheckm(type, checkm.checkm, toolCheckm, toolCheckm, Channel.value(".tsv")) \
	| mix(_wCollectCheckm2(type, checkm2.checkm, toolCheckm2, toolCheckm2, Channel.value(".tsv"))) \
	| set { checkmSelected }

     // Prepare checkm output file
     checkmSelected | splitCsv(sep: '\t', header: true) \
	| map { sample, checkmDict -> checkmDict } \
	| set { checkmList }

     // Prepare gtdb output file
     gtdb.combined | map { sample, gtdb, size -> gtdb } \
	| splitCsv(sep: '\t', header: true) \
	| set { gtdbCombinedList }

     // Prepare missing gtdb output file
     gtdb.missing | map { sample, gtdb, size -> gtdb } \
	| splitCsv(sep: '\t', header: true) \
	| set { gtdbMissingList }

     _wCollectGtdbtkSummary(type, gtdb.summaryRawCombined, Channel.value("gtdb"), \
	Channel.value("gtdbtk"), Channel.value("_summary_raw_combined.tsv")) \
	| set { gtdbCombinedRawSummaryListFiles }

     _wCollectGtdbtkArchea(type, gtdb.archea, Channel.value("gtdb"), \
	Channel.value("gtdbtk"), Channel.value(".ar122.summary.tsv")) \
	| set { gtdbCombinedArchSummaryListFiles }

     _wCollectGtdbtkBac(type, gtdb.bacteria, Channel.value("gtdb"), \
	Channel.value("gtdbtk"), Channel.value(".bac120.summary.tsv")) \
	| set { gtdbCombinedBacSummaryListFiles }

     _wCollectGtdbtkUnclassified(type, gtdb.unclassified, Channel.value("gtdb"), \
	Channel.value("gtdbtk"), Channel.value("_unclassified.tsv")) \
	| set { gtdbCombinedUnclassifiedListFiles }

     _wCollectGtdbtkCombined(type, gtdb.combined, Channel.value("gtdb"), \
	Channel.value("gtdbtk"), Channel.value("_combined.tsv")) \
	| set { gtdbCombinedListFiles }

     _wCollectGtdbtkMissing(type, gtdb.missing, Channel.value("gtdb"), \
	Channel.value("gtdbtk"), Channel.value("_missing_bins.tsv")) \
	| set { gtdbMissingListFiles }

     pGtdbtk.out.logs | mix(pCheckM.out.logs) | mix(pCheckM2.out.logs) | pDumpLogs 

   emit:
     checkm = checkmList
     gtdb =  gtdbCombinedList
     gtdbCombinedSummaryFiles = gtdbCombinedRawSummaryListFiles
     gtdbMissing = gtdbMissingList
     checkmFiles = checkmSelected
}





