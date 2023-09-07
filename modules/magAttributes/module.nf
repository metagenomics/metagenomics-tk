include { wSaveSettingsList } from '../config/module'

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

    publishDir params.output, mode: "${params.publishDirMode}", saveAs: { filename -> getOutput("${sample}", params.runid, "checkm", filename) }, \
      pattern: "{**.tsv}"

    when params.steps.containsKey("magAttributes") && params.steps.magAttributes.containsKey("checkm")

    containerOptions Utils.getDockerMount(params.steps?.magAttributes?.checkm?.database, params)

    beforeScript "mkdir -p ${params.polished.databases}"

    label 'highmemMedium'

    input:
    tuple val(sample), val(ending), path(bins) 

    output:
    tuple path("${sample}_checkm_*.tsv", type: "file"), val("${sample}"), emit: checkm
    tuple env(FILE_ID), val("${output}"), val(params.LOG_LEVELS.INFO), file(".command.sh"), \
	file(".command.out"), file(".command.err"), file(".command.log"), emit: logs

    shell:
    output = getOutput("${sample}", params.runid, "checkm", "")
    S5CMD_PARAMS=params?.steps?.magAttributes?.checkm?.database?.download?.s5cmd?.params ?: "" 
    DOWNLOAD_LINK=params?.steps?.magAttributes?.checkm?.database?.download?.source ?: ""
    MD5SUM=params.steps?.magAttributes?.checkm?.database?.download?.md5sum ?: ""
    EXTRACTED_DB=params.steps?.magAttributes?.checkm?.database?.extractedDBPath ?: ""
    template 'checkm.sh'
    
}


process pGtdbtk {

    container "${params.gtdbtk_image}"

    label 'highmemMedium'

    tag "Sample: $sample"

    publishDir params.output, mode: "${params.publishDirMode}", saveAs: { filename -> getOutput("${sample}",params.runid ,"gtdb", filename) }, \
      pattern: "{**.tsv,**.tree}"

    when params.steps.containsKey("magAttributes") && params.steps.magAttributes.containsKey("gtdb")

    containerOptions Utils.getDockerMount(params?.steps?.magAttributes?.gtdb?.database, params)

    beforeScript "mkdir -p ${params.polished.databases}"

    input:
    tuple val(sample), val(ending), path(bins) 

    output:
    tuple path("chunk_*_${sample}_gtdbtk.bac120.summary.tsv"), val("${sample}"), optional: true, emit: bacteria
    tuple path("chunk_*_${sample}_gtdbtk.ar122.summary.tsv"), val("${sample}"), optional: true, emit: archea
    tuple path("chunk_*_${sample}_gtdbtk_unclassified.tsv"), val("${sample}"), optional: true, emit: unclassified
    tuple path("*.tree"), val("${sample}"), optional: true, emit: tree
    tuple path("chunk_*_${sample}_gtdbtk_combined.tsv"), val("${sample}"), optional: true, emit: combined
    tuple env(FILE_ID), val("${output}"), val(params.LOG_LEVELS.INFO), file(".command.sh"), \
	file(".command.out"), file(".command.err"), file(".command.log"), emit: logs

    shell:
    output = getOutput("${sample}", params.runid, "gtdb", "")
    S5CMD_PARAMS=params?.steps?.magAttributes?.gtdb?.database?.download?.s5cmd?.params ?: ""
    DOWNLOAD_LINK=params?.steps?.magAttributes?.gtdb?.database?.download?.source ?: ""
    MD5SUM=params.steps?.magAttributes?.gtdb?.database?.download?.md5sum ?: ""
    EXTRACTED_DB=params.steps?.magAttributes?.gtdb?.database?.extractedDBPath ?: ""
    GTDB_PARAMS=params.steps.magAttributes.gtdb.additionalParams ?: ""
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
  
     _wMagAttributes(mags)
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
     mags
   main:
     mags | _wMagAttributes
   emit:
     checkm = _wMagAttributes.out.checkm
     gtdb = _wMagAttributes.out.gtdb
}


/*
*
* Method takes a list of the form [SAMPLE, [BIN1 path, BIN2 path]] as input
* and produces a flattend list which is grouped by dataset and sample.
* The output has the form [SAMPLE, file ending (e.g. .fa), [BIN 1 path, BIN 2 path]]
*
*/
def groupBins(binning, buffer){
  def chunkList = [];
  def SAMPLE_IDX = 0;
  def BIN_PATHS_IDX = 1;
  binning[BIN_PATHS_IDX].collate(buffer).each {  
       it.groupBy{ 
            bin -> file(bin).name.substring(file(bin).name.lastIndexOf(".")) 
       }.each {
            ending, group -> chunkList.add([binning[SAMPLE_IDX],  ending, group]);
       }
  }
  return chunkList;
}


workflow _wMagAttributes {
   take: 
     bins
   main:
     GTDB_DEFAULT_BUFFER = 500
     CHECKM_DEFAULT_BUFFER = 30
     BIN_FILES_INPUT_IDX = 1

     DATASET_IDX = 0
     FILE_ENDING_IDX = 1
     BIN_FILES_IDX = 2
     BIN_FILES_OUTPUT_GROUP_IDX = 0
     BIN_FILES_OUTPUT_IDX = 0
     DATASET_OUTPUT_IDX = 1

     // get file ending of bin files (.fa, .fasta, ...) and group by file ending and dataset
     bins | flatMap({n -> groupBins(n, params?.steps?.magAttributes?.checkm?.buffer ?: CHECKM_DEFAULT_BUFFER)}) \
       | pCheckM | set {checkm}
     bins | flatMap({n -> groupBins(n, params?.steps?.magAttributes?.gtdb?.buffer ?: GTDB_DEFAULT_BUFFER )}) \
       | pGtdbtk | set {gtdb}

     // Prepare checkm output file
     checkm.checkm | groupTuple(by: DATASET_OUTPUT_IDX, remainder: true) | map { it -> it[BIN_FILES_OUTPUT_GROUP_IDX] }  | flatten | map { bin -> file(bin) } \
       | collectFile(keepHeader: true, newLine: false ){ item -> [ "bin_attributes.tsv", item.text ] } \
       | splitCsv(sep: '\t', header: true) \
       | set{ checkm_list } 

     if(params.summary){
       // collect checkm files for checkm results across multiple datasets
       checkm.checkm \
          | collectFile(newLine: false, keepHeader: true, storeDir: params.output + "/summary/"){ item ->
         [ "checkm.tsv", item[BIN_FILES_OUTPUT_IDX].text  ]
       }

       // collect gtdb files for results across multiple datasets
       gtdb.bacteria | collectFile(newLine: false, keepHeader: true, storeDir: params.output + "/summary/"){ item ->
         [ "${item[DATASET_IDX]}_bacteria_gtdbtk.tsv", item[BIN_FILES_OUTPUT_IDX].text  ]
       }

       gtdb.archea | collectFile(newLine: false, keepHeader: true, storeDir: params.output + "/summary/"){ item ->
         [ "${item[DATASET_IDX]}_archea_gtdbtk.tsv", item[BIN_FILES_OUTPUT_IDX].text  ]
       }
     }

     pGtdbtk.out.logs | mix(pCheckM.out.logs) | pDumpLogs 

   emit:
     checkm = checkm_list
     gtdb = gtdb.combined
}
