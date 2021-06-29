nextflow.enable.dsl=2

params.ending = ".fa"

MODULE="magAttributes"
VERSION="0.2.0"
def getOutput(SAMPLE, RUNID, TOOL, filename){
    return SAMPLE + '/' + RUNID + '/' + MODULE + '/' + VERSION + '/' + TOOL + '/' + filename
}


process pCmseq {

    container "pbelmann/cmseq:${params.cmseq_tag}"

    label 'tiny'

    errorStrategy 'retry'

    when params.steps.magAttributes.containsKey("prokka")

    input:
    tuple val(sample), file(gff), file(bam), file(bai)

    output:
    file("${sample}_${bin}.txt")

    shell:
    '''
    polymut.py --mincov 10 --gff_file !{gff} !{bam} > !{sample}_!{gff}.txt
    '''
}


process pCheckM {

    container "pbelmann/checkm:${params.checkm_tag}"

    errorStrategy 'ignore'

    publishDir params.output, saveAs: { filename -> getOutput("${sample}", params.runid, "checkm", filename) }

    when params.steps.magAttributes.containsKey("checkm")

    containerOptions " --user 1000:1000  --volume ${params.steps.magAttributes.checkm.database}:/.checkm "

    label 'medium'

    input:
    tuple val(sample), val(ending), path(bins) 

    output:
    tuple path("${sample}_checkm_*.tsv", type: "file"), val("${sample}")

    shell:
    template 'checkm.sh'
    
}


process pGtdbtk {

    container "ecogenomic/gtdbtk:${params.gtdbtk_tag}"

    errorStrategy 'ignore'

    label 'large'

    publishDir params.output, saveAs: { filename -> getOutput("${sample}",params.runid ,"gtdb", filename) }

    when params.steps.magAttributes.containsKey("gtdb")

    containerOptions " --user 1000:1000  --volume ${params.steps.magAttributes.gtdb.database}:/refdata"
   
    input:
    tuple val(sample), val(ending), path(bins) 

    output:
    tuple path("chunk_*_${sample}_gtdbtk.bac120.summary.tsv"), val("${sample}"), optional: true, emit: bacteria
    tuple path("chunk_*_${sample}_gtdbtk.ar122.summary.tsv"), val("${sample}"), optional: true, emit: archea
    tuple path("${sample}_gtdbtk_*.tsv"), val("${sample}"), optional: true, emit: combined

    shell:
    template 'gtdb.sh'
}


process pProkka {

    container "staphb/prokka:${params.prokka_tag}"

    errorStrategy 'ignore'

    label 'small'

    publishDir params.output, saveAs: { filename -> getOutput("${sample}",params.runid ,"prokka", filename) }

    when params.steps.magAttributes.containsKey("prokka")

    input:
    tuple val(sample), file(bin)

    output:
    tuple file("*.gff"), env(BIN_ID), val("${sample}"), emit: gff 
    tuple file("*.err"), env(BIN_ID), val("${sample}"), emit: err 
    tuple file("*.faa"), env(BIN_ID), val("${sample}"), emit: faa 
    tuple file("*.fna"), env(BIN_ID), val("${sample}"), emit: fna 
    tuple file("*.ffn"), env(BIN_ID), val("${sample}"), emit: ffn 
    tuple file("*.fsa"), env(BIN_ID), val("${sample}"), emit: fsa 
    tuple file("*.gbk"), env(BIN_ID), val("${sample}"), emit: gbk
    tuple file("*.log"), env(BIN_ID), val("${sample}"), emit: log
    tuple file("*.tbl"), env(BIN_ID), val("${sample}"), emit: tbl
    tuple file("*.sqn"), env(BIN_ID), val("${sample}"), emit: sqn
    tuple file("*.txt"), env(BIN_ID), val("${sample}"), emit: txt
    tuple file("*.tsv"), env(BIN_ID), val("${sample}"), emit: tsv

    shell:
    '''
    prokka --cpus !{task.cpus} !{bin} --outdir out
    BIN=!{bin}
    BIN_PRAEFIX=$(echo "${BIN%%.*}")
    BIN_ID="${BIN_PRAEFIX}"
    for f in out/* ; do suffix=$(echo "${f#*.}"); mv $f ${BIN_PRAEFIX}.${suffix}; done
    '''
}

def flattenBins(binning){
  def chunkList = [];
  binning[1].each {
     chunkList.add([binning[0], it]);
  }
  return chunkList;
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
*           PATH    DATASET
*           /vol/spool/fragmentRecruitment_20210607/work/3d/9a45b85c15b22a6bb8ed4635391a40/ERR2019981.bam   ERR2019981.bam
*
*/
workflow wCMSeqWorkflowFile {
   take:
      genomes
      alignments
   main:
      wMagAttributesFile(genomes)
      Channel.fromPath(alignments) | splitCsv(sep: '\t', header: true) \
       | map { sample -> [sample.BAM, file(sample.PATH), sample.BAI] } | set {alignments}

      SAMPLE_BINID_IDX = [1,2] 
      SAMPLE_ID_IDX = 1
      GFF_IDX = 2
      ALIGNMENT_IDX = 5
      ALIGNMENT_INDEX_IDX = 5
      wMagAttributesFile.out.prokka_gff \
         | join(wMagAttributesFile.out.prokka_fna, by: SAMPLE_BINID_IDX, remainder: true) \
         | combine(alignments) | map { it -> [it[SAMPLE_ID_IDX],it[GFF_IDX],it[ALIGNMENT_IDX], it[ALIGNMENT_INDEX_IDX] ] } |  pCmseq

}


/**
* This workflow provides the same functionality as wMagAttributesList.
* Input:
* It accepts as input all a tsv table of the following format:
* 
* DATASET	PATH
* SAMPLe_NAME	/path/to/bin.fa
*
*/
workflow wMagAttributesFile {
   take:
     mags
   main:
     DATASET_IDX = 0
     mags | splitCsv(sep: '\t', header: true) \
       | map { sample -> [sample.DATASET, file(sample.PATH)] } \
       | groupTuple(by: DATASET_IDX) | _wMagAttributes
   emit:
     checkm = _wMagAttributes.out.checkm
     prokka_err = _wMagAttributes.out.prokka_err
     prokka_faa = _wMagAttributes.out.prokka_faa
     prokka_ffn = _wMagAttributes.out.prokka_ffn
     prokka_fna = _wMagAttributes.out.prokka_fna
     prokka_fsa = _wMagAttributes.out.prokka_fsa
     prokka_gbk = _wMagAttributes.out.prokka_gbk
     prokka_gff = _wMagAttributes.out.prokka_gff
     prokka_log = _wMagAttributes.out.prokka_log
     prokka_sqn = _wMagAttributes.out.prokka_sqn
     prokka_tbl = _wMagAttributes.out.prokka_tbl
     prokka_tsv = _wMagAttributes.out.prokka_tsv
     prokka_txt = _wMagAttributes.out.prokka_txt


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
*  * Prokka channels for every possible output file err,faa,ffn,fna,fsa,gbk,gff,log,sqn,tbl,tsv,txt,
*    Every prokka channel has the following output format [SAMPLE_NAME, FILE]
*
*/
workflow wMagAttributesList {
   take: 
     mags
   main:
     mags | _wMagAttributes
   emit:
     checkm = _wMagAttributes.out.checkm
     prokka_err = _wMagAttributes.out.prokka_err
     prokka_faa = _wMagAttributes.out.prokka_faa
     prokka_ffn = _wMagAttributes.out.prokka_ffn
     prokka_fna = _wMagAttributes.out.prokka_fna
     prokka_fsa = _wMagAttributes.out.prokka_fsa
     prokka_gbk = _wMagAttributes.out.prokka_gbk
     prokka_gff = _wMagAttributes.out.prokka_gff
     prokka_log = _wMagAttributes.out.prokka_log
     prokka_sqn = _wMagAttributes.out.prokka_sqn
     prokka_tbl = _wMagAttributes.out.prokka_tbl
     prokka_tsv = _wMagAttributes.out.prokka_tsv
     prokka_txt = _wMagAttributes.out.prokka_txt
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
     BIN_FILES_OUTPUT_GROUP_IDX = 1
     BIN_FILES_OUTPUT_IDX = 0

     bins  | flatMap({n -> flattenBins(n)}) | set {binFlattenedList}

     // get file ending of bin files (.fa, .fasta, ...) and group by file ending and dataset
     binFlattenedList | map { it -> def path=file(it[BIN_FILES_INPUT_IDX]); [it[DATASET_IDX], path.name.substring(path.name.lastIndexOf(".")), path]} \
        | set { flattenedListEnding }

     flattenedListEnding  | groupTuple(by: [DATASET_IDX, FILE_ENDING_IDX], size: params?.steps?.magAttributes?.checkm?.buffer ?: CHECKM_DEFAULT_BUFFER, remainder: true) \
        | pCheckM  | set{ checkm }

     flattenedListEnding |  groupTuple(by: [DATASET_IDX, FILE_ENDING_IDX], size: params?.steps?.magAttributes?.gtdb?.buffer ?: GTDB_DEFAULT_BUFFER, remainder: true) \
        | pGtdbtk | set{ gtdb }

     binFlattenedList | pProkka 

     checkm | groupTuple(by: DATASET_IDX, remainder: true) | map { it -> it[BIN_FILES_OUTPUT_GROUP_IDX] }  | flatten | map { bin -> file(bin) } \
       | collectFile(keepHeader: true, newLine: false ){ item -> [ "bin_attributes.tsv", item.text ] } \
       | splitCsv(sep: '\t', header: true) \
       | set{ checkm_list } 

     checkm \
        | collectFile(newLine: false, keepHeader: true, storeDir: params.output + "/summary/"){ item ->
       [ "checkm.tsv", item[BIN_FILES_OUTPUT_IDX].text  ]
     }

     gtdb.bacteria | collectFile(newLine: false, keepHeader: true, storeDir: params.output + "/summary/"){ item ->
       [ "${item[DATASET_IDX]}_bacteria_gtdbtk.tsv", item[BIN_FILES_OUTPUT_IDX].text  ]
     }

     gtdb.archea | collectFile(newLine: false, keepHeader: true, storeDir: params.output + "/summary/"){ item ->
       [ "${item[DATASET_IDX]}_archea_gtdbtk.tsv", item[BIN_FILES_OUTPUT_IDX].text  ]
     }
   emit:
     checkm = checkm_list
     prokka_err = pProkka.out.err
     prokka_faa = pProkka.out.faa
     prokka_ffn = pProkka.out.ffn
     prokka_fna = pProkka.out.fna
     prokka_fsa = pProkka.out.fsa
     prokka_gbk = pProkka.out.gbk
     prokka_gff = pProkka.out.gff
     prokka_log = pProkka.out.log
     prokka_sqn = pProkka.out.sqn
     prokka_tbl = pProkka.out.tbl
     prokka_tsv = pProkka.out.tsv
     prokka_txt = pProkka.out.txt
}
