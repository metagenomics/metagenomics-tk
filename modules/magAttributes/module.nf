nextflow.enable.dsl=2

MODULE="magAttributes"
VERSION="1.0.0"
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

    when params.steps.containsKey("magAttributes") && params.steps.magAttributes.containsKey("checkm")

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

    when params.steps.containsKey("magAttributes") && params.steps.magAttributes.containsKey("gtdb")

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

    container "quay.io/biocontainers/prokka:${params.prokka_tag}"

    errorStrategy 'ignore'

    label 'small'

    time '3h'

    publishDir params.output, saveAs: { filename -> getOutput("${sample}",params.runid ,"prokka", filename) }

    when params.steps.containsKey("magAttributes") && params.steps.magAttributes.containsKey("prokka")

    input:
    tuple val(sample), file(bin), val(domain)

    output:
    tuple file("*.gff.gz"), env(BIN_ID), val("${sample}"), emit: gff 
    tuple file("*.err"), env(BIN_ID), val("${sample}"), emit: err 
    tuple file("*.faa.gz"), env(BIN_ID), val("${sample}"), emit: faa 
    tuple file("*.fna.gz"), env(BIN_ID), val("${sample}"), emit: fna 
    tuple file("*.ffn.gz"), env(BIN_ID), val("${sample}"), emit: ffn 
    tuple file("*.fsa.gz"), env(BIN_ID), val("${sample}"), emit: fsa 
    tuple file("*.gbk.gz"), env(BIN_ID), val("${sample}"), emit: gbk
    tuple file("*.log"), env(BIN_ID), val("${sample}"), emit: log
    tuple file("*.tbl.gz"), env(BIN_ID), val("${sample}"), emit: tbl
    tuple file("*.sqn.gz"), env(BIN_ID), val("${sample}"), emit: sqn
    tuple file("*.txt"), env(BIN_ID), val("${sample}"), emit: txt
    tuple file("*.tsv"), env(BIN_ID), val("${sample}"), emit: tsv

    shell:
    '''
    # Prepare Input Variables
    BIN=!{bin}
    BIN_PREFIX=$(echo "${BIN%.*}")
    BIN_ID="$(basename !{bin})"

    # Run Prokka
    prokka  --cpus !{task.cpus} !{bin}   --outdir out --kingdom !{domain}

    # Prepare output according to magAttributes specification
    for f in out/* ; do suffix=$(echo "${f##*.}"); mv $f ${BIN_PREFIX}.${suffix}; done
    sed -i  -e "2,$ s/^/!{sample}\t${BIN_ID}\t/"  -e "1,1 s/^/SAMPLE\tBIN_ID\t/g" *.tsv
    mv *.tsv !{sample}_prokka_${BIN_ID}.tsv
    gzip --best *gff *.faa *.fna *.ffn *.fsa *.gbk *.sqn *tbl
    '''

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


/*
*
* Method takes a list of the form [SAMPLE, [BIN1 path, BIN2 path]] as input
* and produces a flattend list of the form [SAMPLE, BIN 1 path, BIN 2 path]
*
*/
def flattenBins(binning){
  def chunkList = [];
  def SAMPLE_IDX = 0;
  def BIN_PATHS_IDX = 1;
  binning[BIN_PATHS_IDX].each {
     chunkList.add([binning[SAMPLE_IDX], it]);
  }
  return chunkList;
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
     bins | flatMap({n -> flattenBins(n)}) | set {binFlattenedList}
     bins | flatMap({n -> groupBins(n, params?.steps?.magAttributes?.checkm?.buffer ?: CHECKM_DEFAULT_BUFFER)}) \
       | pCheckM | set {checkm}
     bins | flatMap({n -> groupBins(n, params?.steps?.magAttributes?.gtdb?.buffer ?: GTDB_DEFAULT_BUFFER )}) \
       | pGtdbtk | set {gtdb}

     GTDB_FILE_IDX = 0
     CLASSIFICATION_IDX = 0
     DOMAIN_IDX = 0
     GTDB_OUTPUT_DATASET_IDX = 1
     BIN_FILE_IDX = 1
     BIN_FILE_PROKKA_INPUT_IDX = 2
     DOMAIN_PROKKA_INPUT_IDX = 3

     // if GTDB is enabled then get domain classification for prokka, otherwise just use the user provided default
     if(params?.steps?.magAttributes?.gtdb){

       gtdb.combined | filter(it -> file(it[GTDB_FILE_IDX]).text?.trim()) \
         | splitCsv(sep: '\t', header: true) \
         | map { it ->  def command = it[GTDB_FILE_IDX].classification.split(';')[DOMAIN_IDX].minus('d__'); [it[GTDB_FILE_IDX].SAMPLE, it[GTDB_FILE_IDX].BIN_ID, command] } \
         | set { gtdbDomain }

       binFlattenedList  |  map {it -> [it[DATASET_IDX], file(it[BIN_FILE_IDX]).name, it[BIN_FILE_IDX]]} \
           | join(gtdbDomain, by:[DATASET_IDX, BIN_FILE_IDX], remainder: true) \
           | map { it -> [it[DATASET_IDX], it[BIN_FILE_PROKKA_INPUT_IDX], it[DOMAIN_PROKKA_INPUT_IDX]?:params?.steps?.magAttributes?.prokka?.defaultKingdom  ]} \
           | set { prokkaInput }

     } else {
       binFlattenedList | map { bin  -> [ bin[DATASET_IDX], bin[BIN_FILE_IDX], params?.steps?.magAttributes?.prokka?.defaultKingdom ] } \
           | set { prokkaInput }
     }

     prokkaInput  | pProkka

     // Prepare checkm output file
     checkm | groupTuple(by: DATASET_OUTPUT_IDX, remainder: true) | map { it -> it[BIN_FILES_OUTPUT_GROUP_IDX] }  | flatten | map { bin -> file(bin) } \
       | collectFile(keepHeader: true, newLine: false ){ item -> [ "bin_attributes.tsv", item.text ] } \
       | splitCsv(sep: '\t', header: true) \
       | set{ checkm_list } 

     if(params.summary){
       // collect checkm files for checkm results across multiple datasets
       checkm \
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
