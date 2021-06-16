nextflow.enable.dsl=2

params.ending = ".fa"

MODULE="magAttributes"
VERSION="0.2.0"
def getOutput(SAMPLE, RUNID, TOOL, filename){
    return SAMPLE + '/' + RUNID + '/' + MODULE + '/' + VERSION + '/' + TOOL + '/' + filename
}

process pCheckM {

    container "pbelmann/checkm:${params.checkm_tag}"

    errorStrategy 'ignore'

    publishDir params.output, saveAs: { filename -> getOutput("${sample}", params.runid, "checkm", filename) }

    when params.steps.magAttributes.containsKey("checkm")

    containerOptions " --user 1000:1000  --volume ${params.steps.magAttributes.checkm.database}:/.checkm "

    label 'medium'

    input:
    tuple val(sample), path(bins) 

    output:
    tuple val("${sample}"), path("chunk_*_${sample}_checkm.txt", type: "file") 

    shell:
    '''
    TYPE=checkm_out
    count=`ls -1 *.fasta 2>/dev/null | wc -l`
    if [ $count != 0 ]
    then 
       FILE_ENDING=".fasta"
    else
       FILE_ENDING=".fa"
    fi 

    echo '{"dataRoot": "/.checkm", "remoteManifestURL": "https://data.ace.uq.edu.au/public/CheckM_databases/", "manifestType": "CheckM", "remoteManifestName": ".dmanifest", "localManifestName": ".dmanifest"}' > /tmp/DATA_CONFIG
    mkdir out
    checkm tree -x $FILE_ENDING --reduced_tree --pplacer_threads !{task.cpus}  -t !{task.cpus} -x !{params.ending} . out &> tree.log
    checkm tree_qa out &> tree_qa.log
    checkm lineage_set out out/marker &> lineage.log
    checkm analyze -x $FILE_ENDING -t !{task.cpus} out/marker . out &> analyze.log
    FILE=$(mktemp chunk_XXXXXXXXXX_!{sample}_checkm.txt)
    checkm qa --tab_table -t !{task.cpus} -f checkm.txt out/marker out  &> qa.log

    echo "SAMPLE\tBIN_ID\tMarker lineage\t# genomes\t# markers\t# marker sets\t0\t1\t2\t3\t4\t5+\tCOMPLETENESS\tCONTAMINATION\tHETEROGENEITY" > checkm_tmp.tsv
    tail -n +2 checkm.txt | sed "s/^/!{sample}\t/g"  >> checkm_tmp.tsv

    echo "PATH" > path.tsv
    tail -n +2 checkm.txt | cut -f 1 | sed "s/$/${FILE_ENDING}/g" | xargs -I {} readlink -f {} >> path.tsv

    paste -d$'\t' path.tsv checkm_tmp.tsv > $FILE 
    '''
}


process pGtdbtk {

    container "ecogenomic/gtdbtk:${params.gtdbtk_tag}"

    errorStrategy 'ignore'

    label 'large'

    publishDir params.output, saveAs: { filename -> getOutput("${sample}",params.runid ,"gtdb", filename) }

    when params.steps.magAttributes.containsKey("gtdb")

    containerOptions " --user 1000:1000  --volume ${params.steps.magAttributes.gtdb.database}:/refdata"
   
    input:
    tuple val(sample), path(bins) 

    output:
    tuple val("${sample}"), path("chunk_*_${sample}_gtdbtk.bac120.summary.tsv"), optional: true, emit: bacteria
    tuple val("${sample}"), path("chunk_*_${sample}_gtdbtk.ar122.summary.tsv"), optional: true, emit: archea

    shell:
    '''
    count=`ls -1 *.fasta 2>/dev/null | wc -l`
    if [ $count != 0 ]
    then 
       FILE_ENDING=".fasta"
    else
       FILE_ENDING=".fa"
    fi

    mkdir output
    readlink -f !{bins} > bin.path
    paste -d$'\t' bin.path <(for p in $(cat bin.path); do basename $p; done) > input.tsv
    gtdbtk classify_wf --batchfile input.tsv --out_dir output --cpus !{task.cpus}  --extension ${FILE_ENDING}
    touch output/gtdbtk.bac120.summary.tsv
    touch output/gtdbtk.ar122.summary.tsv
    FILE_BAC=$(mktemp chunk_XXXXXXXXXX_!{sample}_gtdbtk.bac120.summary.tsv)
    FILE_ARC=$(mktemp chunk_XXXXXXXXXX_!{sample}_gtdbtk.ar122.summary.tsv)

    sed "s/^/SAMPLE\t/g" <(head -n 1 output/gtdbtk.bac120.summary.tsv) > $FILE_BAC 
    sed "s/^/!{sample}\t/g"  <(tail -n +2 output/gtdbtk.bac120.summary.tsv) >> $FILE_BAC 

    sed "s/^/SAMPLE\t/g" <(head -n 1 output/gtdbtk.ar122.summary.tsv) > $FILE_ARC 
    sed "s/^/!{sample}\t/g" <(tail -n +2 output/gtdbtk.ar122.summary.tsv) >> $FILE_ARC 
    '''
}


process pProkka {

    container "staphb/prokka:${params.prokka_tag}"

    errorStrategy 'ignore'

    label 'large'

    publishDir params.output, saveAs: { filename -> getOutput("${sample}",params.runid ,"prokka", filename) }

    when params.steps.magAttributes.containsKey("prokka")

    input:
    tuple val(sample), file(bin)

    output:
    tuple file("*.gff"), val("${sample}"), emit: gff 
    tuple file("*.err"), val("${sample}"), emit: err 
    tuple file("*.faa"), val("${sample}"), emit: faa 
    tuple file("*.fna"), val("${sample}"), emit: fna 
    tuple file("*.ffn"), val("${sample}"), emit: ffn 
    tuple file("*.fsa"), val("${sample}"), emit: fsa 
    tuple file("*.gbk"), val("${sample}"), emit: gbk
    tuple file("*.log"), val("${sample}"), emit: log
    tuple file("*.tbl"), val("${sample}"), emit: tbl
    tuple file("*.sqn"), val("${sample}"), emit: sqn
    tuple file("*.txt"), val("${sample}"), emit: txt
    tuple file("*.tsv"), val("${sample}"), emit: tsv

    shell:
    '''
    prokka --cpus !{task.cpus} !{bin} --outdir out
    for f in out/* ; do suffix=$(echo "${f#*.}"); mv $f !{sample}_$(basename !{bin}).${suffix}; done
    '''
}

def flattenBins(binning){
  def chunkList = [];
  binning[1].each {
     chunkList.add([binning[0], it]);
  }
  return chunkList;
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
     DATASET_IDX = 0
     BIN_FILES_IDX = 1

     bins  | flatMap({n -> flattenBins(n)}) | set {binFlattenedList}
     binFlattenedList \
        | groupTuple(by: [DATASET_IDX], size: params?.steps?.magAttributes?.checkm?.buffer ?: CHECKM_DEFAULT_BUFFER, remainder: true) \
        | pCheckM  | set{ checkm }

     binFlattenedList \
        |  groupTuple(by: [DATASET_IDX], size: params?.steps?.magAttributes?.gtdb?.buffer ?: GTDB_DEFAULT_BUFFER, remainder: true) \
        | pGtdbtk | set{ gtdb }

     binFlattenedList | pProkka 

     checkm | groupTuple(by: DATASET_IDX, remainder: true) | map { it -> it[BIN_FILES_IDX] }  | flatten | map { bin -> file(bin) } \
       | collectFile(keepHeader: true, newLine: false ){ item -> [ "bin_attributes.tsv", item.text ] } \
       | splitCsv(sep: '\t', header: true) \
       | set{ checkm_list } 

     checkm \
        | collectFile(newLine: false, keepHeader: true, storeDir: params.output + "/summary/"){ item ->
       [ "checkm.tsv", item[BIN_FILES_IDX].text  ]
     }

     gtdb.bacteria | collectFile(newLine: false, keepHeader: true, storeDir: params.output + "/summary/"){ item ->
       [ "${item[DATASET_IDX]}_bacteria_gtdbtk.tsv", item[BIN_FILES_IDX].text  ]
     }

     gtdb.archea | collectFile(newLine: false, keepHeader: true, storeDir: params.output + "/summary/"){ item ->
       [ "${item[DATASET_IDX]}_archea_gtdbtk.tsv", item[BIN_FILES_IDX].text  ]
     }
     checkm_list | view()
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
