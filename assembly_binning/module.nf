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
params.checkm_database = ""
params.gtdb = false
params.skip = false
params.buffer = 30
params.ending = ".fa"

process runMetaSpades {

    label 'large'

    tag "$sample"

    container 'quay.io/biocontainers/spades:${metaspades_tag}'

    publishDir "${params.output}/${sample}/assembly/metaspades/${metaspades_tag}" 

    when params.metaspades

    input:
    tuple val(sample), path(fastqs, stageAs: 'fastq.fq.gz')

    output:
    tuple val("${sample}"), env(TYPE), path("${sample}/contigs.fasta"), emit: contigs

    shell:
    '''
    metaspades.py --threads 28  --12 fastq.fq.gz -o !{sample}
    TYPE="metaspades"
    '''
}


process runMegahit {

    label 'large'

    tag "$sample"

    container "vout/megahit:${megahit_tag}"

    publishDir "${params.output}/${sample}/assembly/megahit/${megahit_tag}" 

    when params.containsKey("megahit") && params.assembly == "megahit"

    cpus 28
    input:
    tuple val(sample), path(fastqs, stageAs: 'fastq.fq.gz')

    output:
    tuple val("${sample}"), env(TYPE), path("${sample}/final.contigs.fa"), emit: contigs

    shell:
    '''
    megahit -t !{task.cpus} --12 !{fastqs}
    TYPE="megahit" 
    mkdir !{sample}
    mv megahit_out/final.contigs.fa !{sample}/
    '''
}



process runBowtie {

    container "pbelmann/bowtie2:${params.bowtie_tag}"

    label 'large'

    tag "$sample"

    publishDir "${params.output}/${sample}/mapping/bowtie/${params.bowtie_tag}" 

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

    container "metabat/metabat:${params.metabat_tag}"

    errorStrategy 'ignore'

    tag "$sample"

    label 'large'

    publishDir "${params.output}/${sample}/binning/metabat/${params.metabat_tag}" 

    when params.containsKey("binning") && params.binning == "metabat"

    input:
    tuple val(sample), val(TYPE), path(contigs), path(bam)

    output:
    tuple val("${sample}"), env(NEW_TYPE), file("${TYPE}/bin*"), optional: true, emit: bins
    tuple val("${sample}"), val("${TYPE}"), file("${TYPE}/bin*"), optional: true, emit: bins_assembler

    shell:
    '''
    NEW_TYPE="!{TYPE}_metabat"
    runMetaBat.sh !{contigs} !{bam}
    mkdir !{TYPE}
    mv $(basename !{contigs})*/bin* !{TYPE}
    '''
}


process runMaxBin {

    container "quay.io/biocontainers/maxbin2:${params.maxbin_tag}"

  //  errorStrategy 'ignore'

    label 'large'

    tag "$sample"

    when params.maxbin

    publishDir "${params.output}/${sample}/binning/maxbin/${params.maxbin_tag}" 

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

    container "pbelmann/checkm:${params.checkm_tag}"

    errorStrategy 'ignore'

    publishDir "${params.output}/${sample}/checkm/${params.checkm_tag}" 

    when params.postprocessing.containsKey("checkm")

    containerOptions " --user 1000:1000  --volume ${params.postprocessing.checkm.database}:/.checkm "

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
    tail -n +2 checkm.txt | sed "s/^/!{sample}\t!{sample}_/g"  >> checkm_tmp.tsv

    echo "PATH" > path.tsv
    tail -n +2 checkm.txt | cut -f 1 | sed "s/$/${FILE_ENDING}/g" | xargs -I {} readlink -f {} >> path.tsv

    paste -d$'\t' path.tsv checkm_tmp.tsv > $FILE 
    '''
}


process runGtdbtk {

    container "ecogenomic/gtdbtk:${params.gtdbtk_tag}"

    errorStrategy 'ignore'

    label 'medium'

    scratch false

    publishDir "${params.output}/${sample}/gtdb/${params.gtdbtk_tag}" 

    when params.postprocessing.containsKey("gtdb")

    containerOptions " --user 1000:1000  --volume ${params.postprocessing.gtdb.database}:/refdata"
   
    input:
    tuple val(sample), val(TYPE), path(bins, stageAs: "bin*.fa") 

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
    ls -1 bin*.fa > bin.id
    readlink -f bin*.fa > bin.path
    paste -d$'\t' bin.path bin.id > input.tsv
    gtdbtk classify_wf --batchfile input.tsv --out_dir output --cpus !{task.cpus}  --extension ${FILE_ENDING}
    touch output/gtdbtk.bac120.summary.tsv
    touch output/gtdbtk.ar122.summary.tsv
    FILE_BAC=$(mktemp chunk_XXXXXXXXXX_!{sample}_!{TYPE}_gtdbtk.bac120.summary.tsv)
    FILE_ARC=$(mktemp chunk_XXXXXXXXXX_!{sample}_!{TYPE}_gtdbtk.ar122.summary.tsv)
    sed "s/^/!{sample}\t/g" output/gtdbtk.bac120.summary.tsv > $FILE_BAC 
    sed "s/^/!{sample}\t/g" output/gtdbtk.ar122.summary.tsv > $FILE_ARC 
    '''
}




process prokka {

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

process getMappingQuality {

    container 'quay.io/biocontainers/samtools:1.12--h9aed4be_1'

    tag "$sample"

    publishDir "${params.output}/${sample}/mapping/bowtie/"

    errorStrategy 'retry'

    label 'tiny'

    input:
    tuple val(sample), val(TYPE), path(bam) 

    output:
    tuple val("${sample}"), val("${TYPE}"), file("${sample}_${TYPE}_flagstat.tsv"), emit: flagstat_raw
    tuple val("${sample}"), val("${TYPE}"), file("${sample}_${TYPE}_flagstat_passed.tsv"), emit: flagstat_passed
    tuple val("${sample}"), val("${TYPE}"), file("${sample}_${TYPE}_flagstat_failed.tsv"), emit: flagstat_failed

    shell:
    template 'mapping_quality.sh'
}


process runGetReads {

    //container 'biocontainers/samtools:v1.7.0_cv4'

    publishDir "${params.output}/${sample}/reads/bins/" 

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


process runTrimmomatic {

    label 'large'

    tag "$sample"

    publishDir "${params.output}/${sample}/reads/" 

    container "quay.io/biocontainers/trimmomatic:${params.trimmomatic_tag}"

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



process runFastp {

    label 'large'

    tag "$sample"

    publishDir "${params.output}/${sample}/reads/fastp/${params.fastp_tag}"

    container "quay.io/biocontainers/fastp:${params.fastp_tag}"

    input:
    tuple val(sample), path(genomeReads1), path(genomeReads2)

    output:
    tuple val("${sample}"), path("*_R1.p.fastq.gz"), path("*_R2.p.fastq.gz"), emit: fastq
    path("*_fastp_report.html"), emit: html_report
    path("*_fastp_report.json"), emit: json_report
    path("fastp_summary"), emit: fastp_summary
    

    script:
    fq_1_paired = sample + '_R1.p.fastq.gz'
    fq_1_unpaired = sample + '_R1.s.fastq.gz'
    fq_2_paired = sample + '_R2.p.fastq.gz'
    fq_2_unpaired = sample + '_R2.s.fastq.gz'
    template 'fastp.sh'
}


def bufferMetabatSamtools(metabat){
  def chunkList = [];
  metabat[2].each {
     chunkList.add([metabat[0],metabat[1], metabat[3], it]);
  }
  return chunkList;
}

def bufferBins(binning){
  def chunkList = [];
  binning[2].each {
     chunkList.add([binning[0],binning[1], it]);
  }
  return chunkList;
}

process runBBMapInterleave {

    label 'large'

    tag "$sample"

    publishDir "${params.output}/${sample}/reads/bbmap/" 

    container "quay.io/biocontainers/bbmap:${params.bbmap_tag}"

    input:
    tuple val(sample), path(read1, stageAs: "read1.fq.gz"), path(read2, stageAs: "read2.fq.gz")

    output:
    tuple val("${sample}"), path("reads.fq.gz")

    shell:
    """
    reformat.sh in1=read1.fq.gz in2=read2.fq.gz out=reads.fq.gz
    """
}

process runMegahitInterleaved {

    label 'large'

    tag "$sample"

    publishDir "${params.output}/${sample}/assembly/megahit/${params.megahit_tag}" 

    when params.containsKey("megahit") && params.assembly == "megahit"

    errorStrategy 'ignore'

    container "vout/megahit:release-v1.2.9:${params.megahit_tag}"

    input:
    tuple val(sample), path(interleaved_reads, stageAs: "interleaved.fq.gz")

    output:
    tuple val("${sample}"), env(TYPE), path("final.contigs.fa"), emit: contigs
    tuple val("${sample}"), path("interleaved.fastp.fq.gz"), emit: reads_processed
    tuple val("${sample}"), path("fastp_summary.tsv"), emit: fastp_summary


    shell:
    template 'megahit_interleaved.sh'

}


process runMegahitSplit {

    label 'large'

    tag "$sample"

    publishDir "${params.output}/${sample}/assembly/megahit/${params.megahit_tag}" 

    errorStrategy 'ignore'

    when params.containsKey("megahit") && params.assembly == "megahit"

    container "vout/megahit:${params.megahit_tag}"

    input:
    tuple val(sample), path(read1, stageAs: "read1.fq.gz"), path(read2, stageAs: "read2.fq.gz")

    output:
    tuple val("${sample}"), env(TYPE), path("final.contigs.fa"), emit: contigs
    tuple val("${sample}"), path("interleaved.fastp.fq.gz"), emit: reads_processed
    tuple val("${sample}"), path("fastp_summary.tsv"), emit: fastp_summary

    shell:
    template 'megahit_split.sh'

}

process runBBMapDeinterleave {

    label 'large'

    tag "$sample"

    publishDir "${params.output}/${sample}/reads/bbmap/${params.bbmap_tag}" 

    container "quay.io/biocontainers/bbmap:${params.bbmap_tag}"

    input:
    tuple val(sample), path(interleaved_reads, stageAs: "interleaved.fq.gz")

    output:
    tuple val("${sample}"), path("read1.fq.gz"), path("read2.fq.gz")

    shell:
    """
    reformat.sh in=interleaved.fq.gz out1=read1.fq.gz out2=read2.fq.gz
    """
}


workflow binning {
   take: 
     contigs
     input_reads
   main:
     contigs | join(input_reads | mix(input_reads), by: 0) | runBowtie | set { bam }
     bam | getMappingQuality 
     
     getMappingQuality.out.flagstat_passed | collectFile(newLine: false, keepHeader: true, storeDir: params.output ){ item ->
       [ "flagstat_passed.tsv", item[2].text  ]
     }

     getMappingQuality.out.flagstat_failed | collectFile(newLine: false, keepHeader: true, storeDir: params.output ){ item ->
       [ "flagstat_failed.tsv", item[2].text  ]
     }

     contigs | join(bam, by: [0, 1]) | runMetabat | set { metabat }
     contigs | join(input_reads | mix(input_reads), by: 0) | runMaxBin | set { maxbin }

    // metabat.bins | mix(maxbin.bins) | set {bins}
     postprocess(metabat.bins, bam)
   emit:
     postprocess.out

}


workflow postprocess {
   take: 
     bins
     bam
   main:
     bins | flatMap({n -> bufferBins(n)}) | groupTuple(by: [0,1], size: params.buffer, remainder: true) | set{bufferedBins}

     bufferedBins | runCheckM  | set{ checkm }
     bufferedBins | runGtdbtk | set{ gtdb}

     checkm | collectFile(newLine: false, keepHeader: true, storeDir: params.output ){ item ->
       [ "${item[1]}_checkm.tsv", item[0].text  ]
     }

     checkm  | groupTuple(by: 2, remainder: true) | map { it -> it[0] }  | flatten | map { it -> file(it) } | collectFile(keepHeader: true, newLine: false ){ item -> [ "bin_attributes.tsv", item.text ] } | splitCsv(sep: '\t', header: true)  \
       | map { it -> { it.put('COVERAGE', 0); it} } | map { it -> { it.put('N50', 0); it} }  | set{ bins_info } 

     gtdb.bacteria | collectFile(newLine: false, keepHeader: true, storeDir: params.output ){ item ->
       [ "${item[1]}_bacteria_gtdbtk.tsv", item[0].text  ]
     } 

     gtdb.archea | collectFile(newLine: false, keepHeader: true, storeDir: params.output ){ item ->
       [ "${item[1]}_archea_gtdbtk.tsv", item[0].text  ]
     }

     bins | join(bam, by: [0,1]) |  flatMap({n -> bufferMetabatSamtools(n)})  | runGetReads

   emit:
     bins_info

}


workflow assembly_binning {
   take: 
     input_split_raw_reads
   main:
     input_split_raw_reads | runFastp 
     runFastp.out.fastq | runBBMapInterleave | set{input_reads}

     runFastp.out.fastp_summary | collectFile(newLine: false, keepHeader: true, storeDir: params.output ){ item ->
          [ "fastp_summary.tsv", item[1].text  ]
     }
      
     input_reads | runMegahit | set { megahit }
     input_reads | runMetaSpades | set { metaspades }

     metaspades.contigs | mix(megahit.contigs) |  set{contigs}
     binning(contigs, input_reads)
   emit:
     binning.out
}


workflow assembly_binning_input {
     take:
       reads
     main:
       if(params.mode.interleaved){
          if(params.mode.skip){
            reads | splitCsv(sep: '\t', header: true) \
             | map { it -> [ it.SAMPLE, it.READS ]} \
             | runMegahitInterleaved 
            binning(runMegahitInterleaved.out.contigs, runMegahitInterleaved.out.reads_processed)

            runMegahitInterleaved.out.fastp_summary | collectFile(newLine: false, keepHeader: true, storeDir: params.output ){ item ->
              [ "fastp_summary.tsv", item[1].text  ]
            }

          } else {
            reads | splitCsv(sep: '\t', header: true) \
             | map { it -> [ it.SAMPLE, it.READS ]} \
             | runBBMapDeinterleave | assembly_binning
          }
       }

       if(!params.mode.interleaved){
          if(params.mode.skip){
            reads | splitCsv(sep: '\t', header: true) \
             | map { it -> [ it.SAMPLE, it.READS1, it.READS2 ]} \
             | runMegahitSplit 
             
            runMegahitSplit.out.reads_processed | set { processed_reads}
            binning(runMegahitSplit.out.contigs, runMegahitSplit.out.reads_processed)
            runMegahitSplit.out.fastp_summary | collectFile(newLine: false, keepHeader: true, storeDir: params.output ){ item ->
              [ "fastp_summary.tsv", item[1].text ]
            }
          } else {
            reads | splitCsv(sep: '\t', header: true) \
             | map { it -> [ it.SAMPLE, it.READS1, it.READS2 ]} \
             | assembly_binning
          }
       }
    emit:
      bins = binning.out
      processed_reads = processed_reads
}


workflow assembly_binning_input_sra  {
     if(params.sra){
        Channel.fromSRA(['SRR6820513'], apiKey: '9b9acc33c35f7283d76c63eb407a849d1608')  | set{ input_reads }
     }
     assembly_binning(input_reads)
}
