nextflow.enable.dsl=2

params.ending = ".fa"

process pCheckM {

    container "pbelmann/checkm:${params.checkm_tag}"

    errorStrategy 'ignore'

    publishDir "${params.output}/${sample}/checkm/${params.checkm_tag}" 

    when params.steps.magAttributes.containsKey("checkm")

    containerOptions " --user 1000:1000  --volume ${params.steps.magAttributes.checkm.database}:/.checkm "

    label 'medium'

    input:
    tuple val(sample), val(TYPE), path(bins) 

    output:
    tuple path("chunk_*_${sample}_${TYPE}_checkm.txt", type: "file"), val("${sample}"),  val("${TYPE}")

    shell:
    '''
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
    FILE=$(mktemp chunk_XXXXXXXXXX_!{sample}_!{TYPE}_checkm.txt)
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

    publishDir "${params.output}/${sample}/gtdb/${params.gtdbtk_tag}" 

    when params.steps.magAttributes.containsKey("gtdb")

    containerOptions " --user 1000:1000  --volume ${params.steps.magAttributes.gtdb.database}:/refdata"
   
    input:
    tuple val(sample), val(TYPE), path(bins) 

    output:
    tuple path("chunk_*_${sample}_${TYPE}_gtdbtk.bac120.summary.tsv"), val("${TYPE}"), optional: true, emit: bacteria
    tuple path("chunk_*_${sample}_${TYPE}_gtdbtk.ar122.summary.tsv"), val("${TYPE}"), optional: true, emit: archea

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
    FILE_BAC=$(mktemp chunk_XXXXXXXXXX_!{sample}_!{TYPE}_gtdbtk.bac120.summary.tsv)
    FILE_ARC=$(mktemp chunk_XXXXXXXXXX_!{sample}_!{TYPE}_gtdbtk.ar122.summary.tsv)

    sed "s/^/SAMPLE\t/g" <(head -n 1 output/gtdbtk.bac120.summary.tsv) > $FILE_BAC 
    sed "s/^/!{sample}\t/g"  <(tail -n +2 output/gtdbtk.bac120.summary.tsv) >> $FILE_BAC 

    sed "s/^/SAMPLE\t/g" <(head -n 1 output/gtdbtk.ar122.summary.tsv) > $FILE_ARC 
    sed "s/^/!{sample}\t/g" <(tail -n +2 output/gtdbtk.ar122.summary.tsv) >> $FILE_ARC 
    '''
}


process pProkka {

    container "staphb/prokka:${prokka_tag}"

//    publishDir "${params.output}/${sample}/prokka/${prokka_tag}" 

//    errorStrategy 'retry'

    label 'large'

    input:
    tuple file(bam), file(bai), file(binFile), sampleName 

    output:
    tuple file("out/*.gff"), file("${binFile}"), file("${bam}"), file("${bai}"), val("${sampleName}") 

    shell:
    '''
    prokka --cpus 0 !{binFile} --outdir out
    '''
}


process pGetReads {

    //container 'biocontainers/samtools:v1.7.0_cv4'

   // publishDir "${params.output}/${sample}/reads/bins/" 

    tag "$sample"

    errorStrategy 'retry'

    when params.getreads

    label 'tiny'

    input:
    tuple val(sample), val(TYPE), path(bam), path(bin) 

    output:
    file("${sample}_${TYPE}_${bin}.fq.gz")
    file("${sample}_${TYPE}_${bin}.fq.depth")
    file("${sample}_${TYPE}_${bin}.fq.depth.coverage")

    shell:
    '''
    CONTIGS=$(grep ">" !{bin} | tr -d '>' | tr '\n' ' ')
    DEPTH_FILE=!{sample}_!{TYPE}_!{bin}.fq.depth
    samtools view -b !{bam} | samtools sort -o !{bam}.sorted
    samtools index -@ 1 !{bam}.sorted
    samtools view -h -f 3 -b !{bam}.sorted $CONTIGS | samtools sort -n | samtools bam2fq -N -  | gzip --best > !{sample}_!{TYPE}_!{bin}.fq.gz
    samtools view    -f 3 -b !{bam}.sorted $CONTIGS | samtools depth -a - > !{sample}_!{TYPE}_!{bin}.fq.depth
    ALL_POS=$(cut -f 3 $DEPTH_FILE  | wc -l)
    DEPTH=$(cut -f 3 $DEPTH_FILE | paste -sd+ | head | bc)
    COVERAGE_DEPTH=$(bc <<< "scale = 5; $DEPTH / $ALL_POS")
    echo $COVERAGE_DEPTH > !{sample}_!{TYPE}_!{bin}.fq.depth.coverage
    '''
}


def bufferMetabatSamtools(metabat){
  def chunkList = [];
  metabat[2].each {
     chunkList.add([metabat[0],metabat[1], metabat[3], it]);
  }
  return chunkList;
}


def flattenBins(binning){
  def chunkList = [];
  binning[2].each {
     chunkList.add([binning[0],binning[1], it]);
  }
  return chunkList;
}


workflow wMagAttributes {
   take: 
     bins
     bam
   main:
     bins | flatMap({n -> flattenBins(n)}) | set {binList}
     binList | groupTuple(by: [0,1], size: params.steps.magAttributes.checkm.buffer, remainder: true) \
        | pCheckM  | set{ checkm }

     binList |  groupTuple(by: [0,1], size: params.steps.magAttributes.gtdb.buffer, remainder: true) \
        | pGtdbtk | set{ gtdb }

     checkm | collectFile(newLine: false, keepHeader: true, storeDir: params.output + "/summary/"){ item ->
       [ "checkm.tsv", item[0].text  ]
     }

     checkm | groupTuple(by: 2, remainder: true) | map { it -> it[0] }  | flatten | map { it -> file(it) } \
       | collectFile(keepHeader: true, newLine: false ){ item -> [ "bin_attributes.tsv", item.text ] } \
       | splitCsv(sep: '\t', header: true) | set{ bins_info } 

     gtdb.bacteria | collectFile(newLine: false, keepHeader: true, storeDir: params.output + "/summary/"){ item ->
       [ "${item[1]}_bacteria_gtdbtk.tsv", item[0].text  ]
     } 

     gtdb.archea | collectFile(newLine: false, keepHeader: true, storeDir: params.output + "/summary/"){ item ->
       [ "${item[1]}_archea_gtdbtk.tsv", item[0].text  ]
     }

    // bins | join(bam, by: [0,1]) |  flatMap({n -> bufferMetabatSamtools(n)})  | pGetReads

   emit:
     bins_info

}
