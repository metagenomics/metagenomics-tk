nextflow.enable.dsl=2

params.output = "out"
params.input = "sample_all.tsv"
params.quast = false
params.spades = false
params.metaspades = false
params.megahit = false
params.checkm = false
params.metabat = false
params.maxbin = false
params.getreads = false
params.sra = false
params.interleaved = false
params.deinterleaved = false
params.gtdb_database = ""
params.gtdb = false


process runMetaSpades {

    container 'quay.io/biocontainers/spades:3.15.2--h95f258a_1'

    scratch "/vol/scratch"

    publishDir "${params.output}/${sample}/assembly/metaspades/" 

    cpus 28

    when params.metaspades

    input:
    tuple val(sample), path(fastqs, stageAs: 'fastq.fq.gz')

    output:
    tuple val("${sample}"), env(TYPE), path("${sample}/contigs.fasta"), emit: metabat

    shell:
    '''
    metaspades.py --threads 28  --12 fastq.fq.gz -o !{sample}
    TYPE="metaspades"
    '''
}


process runMegahit {

    container 'vout/megahit:release-v1.2.9'

    scratch "/vol/scratch"

    publishDir "${params.output}/${sample}/assembly/megahit/" 

    when params.megahit

    cpus 28
    input:
    tuple val(sample), path(fastqs, stageAs: 'fastq.fq.gz')

    output:
    tuple val("${sample}"), env(TYPE), path("${sample}/final.contigs.fa"), emit: metabat

    shell:
    '''
    megahit -t !{task.cpus} --12 !{fastqs}
    TYPE="megahit" 
    mkdir !{sample}
    mv megahit_out/final.contigs.fa !{sample}/
    '''
}



process runBowtie {

    container 'pbelmann/bowtie2:0.11.0'

    maxForks 40

    scratch "/vol/scratch"

    publishDir "${params.output}/${sample}/mapping/bowtie/" 

    errorStrategy 'retry'

    cpus 28

    input:
    tuple val(sample), val(TYPE), path(contigs), path(fastqs, stageAs: 'fastq.fq.gz')

    output:
    tuple val("${sample}"), val("${TYPE}"), file("${TYPE}_${sample}.bam"), emit: bam

    shell:
    '''
    INDEX=!{sample}.index
    bowtie2-build --threads 28 --quiet !{contigs} $INDEX 
    bowtie2 -p 28 --quiet --very-sensitive -x $INDEX --interleaved fastq.fq.gz | samtools view -F 3584 --threads 28 -bS - | samtools sort --threads 28 - > !{TYPE}_!{sample}.bam
    '''
}


process runMetabat {

    container 'metabat/metabat:v2.15-4-ga101cde'

    scratch "/vol/scratch"

    errorStrategy 'ignore'

    publishDir "${params.output}/${sample}/binning/metabat/" 

    when params.metabat

    cpus 28

    input:
    tuple val(sample), val(TYPE), path(contigs), path(bam)

    output:
    tuple val("${sample}"), env(NEW_TYPE), file("${TYPE}_${sample}/bin*"), optional: true, emit: bins
    tuple val("${sample}"), val("${TYPE}"), file("${TYPE}_${sample}/bin*"), optional: true, emit: bins_assembler

    shell:
    '''
    NEW_TYPE="!{TYPE}_metabat"
    runMetaBat.sh !{contigs} !{bam}
    mkdir !{TYPE}_!{sample}
    mv $(basename !{contigs})*/bin* !{TYPE}_!{sample}
    '''
}


process runMaxBin {

    container 'quay.io/biocontainers/maxbin2:2.2.7--he1b5a44_1'

    scratch "/vol/scratch"
  //  errorStrategy 'ignore'

    cpus 28

    when params.maxbin

    publishDir "${params.output}/${sample}/binning/maxbin/" 

    input:
    tuple val(sample), val(TYPE), path(contigs), path(reads)

    output:
    tuple val("${sample}"), env(NEW_TYPE), file("${TYPE}_${sample}/out.*.fasta"), optional: true, emit: bins

    shell:
    '''
    NEW_TYPE="!{TYPE}_maxbin"
    mkdir !{TYPE}_!{sample}
    run_MaxBin.pl -preserve_intermediate -contig !{contigs} -reads !{reads} -thread 28 -out !{TYPE}_!{sample}/out
    '''
}


process runCheckM {

    container 'pbelmann/checkm:0.12.0'

    errorStrategy 'retry'

    publishDir "${params.output}/${sample}/checkm/" 

    when params.checkm

    containerOptions " --user 1000:1000  --volume ${params.database}:/.checkm "
   
    scratch "/vol/scratch"

    cpus 14

    input:
    tuple val(sample), val(TYPE), path(bins) 

    output:
    tuple path("${sample}_${TYPE}_checkm.txt", type: "file"), val("${TYPE}")

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
    DIRECTORY=$(mktemp -d --suffix=.checkm_out  -p .)   
    checkm qa --tab_table -t !{task.cpus} -f ${DIRECTORY}/!{sample}_!{TYPE}_checkm.txt out/marker out  &> qa.log
    sed  "s/^/!{sample}\t/g" ${DIRECTORY}/!{sample}_!{TYPE}_checkm.txt > !{sample}_!{TYPE}_checkm.txt
    '''
}


process runGtdbtk {

    container 'ecogenomic/gtdbtk:1.4.1'

    errorStrategy 'ignore'

    publishDir "${params.output}/${sample}/gtdb/" 

    when params.gtdb

    containerOptions " --user 1000:1000  --volume ${params.gtdb_database}:/refdata"
   
    scratch "/vol/scratch"

    cpus 28

    input:
    tuple val(sample), val(TYPE), path(bins, stageAs: "bin*.fa") 

    output:
    tuple path("${sample}_${TYPE}_gtdbtk.bac120.summary.tsv"), val("${TYPE}")

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
     
    ls -1 bin*.fa > bin.id
    readlink -f bin*.fa > bin.path
    paste -d$'\t' bin.path bin.id > input.tsv
    gtdbtk classify_wf --batchfile input.tsv --out_dir output --cpus !{task.cpus}  --extension ${FILE_ENDING}
    sed  "s/^/!{sample}\t/g" output/gtdbtk.bac120.summary.tsv > !{sample}_!{TYPE}_gtdbtk.bac120.summary.tsv 
    '''
}



process prokka {

    container 'staphb/prokka:1.14.5'

    scratch "/vol/scratch"

//    publishDir "${params.output}/${sample}/prokka/" 

//    errorStrategy 'retry'

//    time '3h'

    cpus 28

    input:
    tuple file(bam), file(bai), file(binFile), sampleName 

    output:
    tuple file("out/*.gff"), file("${binFile}"), file("${bam}"), file("${bai}"), val("${sampleName}") 

    shell:
    '''
    prokka --cpus 0 !{binFile} --outdir out
    '''
}



process runGetReads {

    //container 'biocontainers/samtools:v1.7.0_cv4'

    publishDir "${params.output}/${sample}/reads/bins/" 

    errorStrategy 'retry'

    scratch "/vol/scratch"

    when params.getreads

    cpus 1

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


process runTrimmomatic {

    cpus 28

    publishDir "${params.output}/${sample}/reads/" 

    container 'quay.io/biocontainers/trimmomatic:0.39--hdfd78af_2'

//    scratch "/vol/scratch"

    input:
    tuple val(sample), path(genomeReads1), path(genomeReads2)

    output:
    tuple val("${sample}"), path("*_R1.p.fastq.gz"), path("*_R2.p.fastq.gz")

    script:

    fq_1_paired = sample + '_R1.p.fastq.gz'
    fq_1_unpaired = sample + '_R1.s.fastq.gz'
    fq_2_paired = sample + '_R2.p.fastq.gz'
    fq_2_unpaired = sample + '_R2.s.fastq.gz'

    """
    trimmomatic \
    PE -phred33 \
    -threads ${task.cpus} \
    ${genomeReads1} \
    ${genomeReads2} \
    $fq_1_paired \
    $fq_1_unpaired \
    $fq_2_paired \
    $fq_2_unpaired \
    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36
    """
}


def bufferMetabatSamtools(metabat){
  def chunkList = [];
  metabat[2].each {
     chunkList.add([metabat[0],metabat[1], metabat[3], it]);
  }
  return chunkList;
}

process runBBMapInterleave {

    cpus 28

    publishDir "${params.output}/${sample}/reads/" 

//    scratch "/vol/scratch"

    container 'quay.io/biocontainers/bbmap:38.90--he522d1c_1'

    input:
    tuple val(sample), path(read1, stageAs: "read1.fq.gz"), path(read2, stageAs: "read2.fq.gz")

    output:
    tuple val("${sample}"), path("reads.fq.gz")

    shell:
    """
    reformat.sh in1=read1.fq.gz in2=read2.fq.gz out=reads.fq.gz
    """
}

process runBBMapDeinterleave {

    cpus 28

    publishDir "${params.output}/${sample}/reads/" 

//    scratch "/vol/scratch"

    container 'quay.io/biocontainers/bbmap:38.90--he522d1c_1'

    input:
    tuple val(sample), path(interleaved_reads, stageAs: "interleaved.fq.gz")

    output:
    tuple val("${sample}"), path("read1.fq.gz"), path("read2.fq.gz")

    shell:
    """
    reformat.sh in=interleaved.fq.gz out1=read1.fq.gz out2=read2.fq.gz
    """
}


workflow assembly_binning {
   take: 
     input_split_raw_reads
   main:
     input_split_raw_reads | runTrimmomatic | runBBMapInterleave | set{input_reads}
      
     input_reads | runMegahit | set { megahit }
     input_reads | runMetaSpades | set { metaspades }

     metaspades.metabat | mix(megahit.metabat) |  set{contigs}
     contigs | join(input_reads | mix(input_reads), by: 0) | runBowtie | set { bam }

     contigs | join(bam, by: [0, 1]) | runMetabat | set { metabat }
     contigs | join(input_reads | mix(input_reads), by: 0) | runMaxBin | set { maxbin }

     metabat.bins | mix(maxbin.bins) | set {bins}

     bins | runCheckM | set{ checkm }
     bins | runGtdbtk

//     checkm | collectFile(skip: 1, newLine: false, keepHeader: true, storeDir: params.output + "/checkm/"){ item ->
//       [ "${item[1]}.txt", item[0].text + '\n' ]
//     }

     metabat.bins_assembler | join(bam, by: [0,1]) |  flatMap({n -> bufferMetabatSamtools(n)})  | runGetReads
}

workflow assembly_binning_input {
     take:
       reads
     main:
       if(params.interleaved){
          reads | splitCsv(sep: '\t', header: true) \
           | map { it -> [ it.SAMPLE, it.READS ]} \
           | runBBMapDeinterleave | set{ input_reads }
       }

       if(params.deinterleaved){
          reads | splitCsv(sep: '\t', header: true) \
           | map { it -> [ it.SAMPLE, it.READS1, it.READS2 ]} \
           | set{ input_reads }
       }

       assembly_binning(input_reads)
}

workflow assembly_binning_input_sra  {
     if(params.sra){
        Channel.fromSRA(['SRR6820513'], apiKey: '9b9acc33c35f7283d76c63eb407a849d1608') | view() | set{ input_reads }
     }
     assembly_binning(input_reads)
}
