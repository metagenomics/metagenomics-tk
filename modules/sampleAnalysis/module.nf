nextflow.enable.dsl=2

def flattenBins(binning){
  def chunkList = [];
  binning[1].each {
     chunkList.add([binning[0], it]);
  }
  return chunkList;
}


MODULE="sampleAnalysis"
VERSION="0.1.0"
def getOutput(SAMPLE, RUNID, TOOL, filename){
    return SAMPLE + '/' + RUNID + '/' + MODULE + '/' + VERSION + '/' + TOOL + '/' + filename
}

process pBowtie {

    container "pbelmann/bowtie2:${params.bowtie_tag}"

    label 'large'

    tag "$sample"

    publishDir params.output, saveAs: { filename -> getOutput("${sample}", params.runid , "bowtie", filename) }

    when params?.steps?.sampleAnalysis?.bowtie != null

    errorStrategy 'retry'

    input:
    tuple val(sample), val(TYPE), path(contigs), path(fastqs, stageAs: 'fastq.fq.gz')

    output:
    tuple val("${sample}"), val("${TYPE}"), file("interleaved.discordand.fq.gz"), optional: true, emit: unmappedReads
    tuple val("${sample}"), val("${TYPE}"), file("${sample}_bowtie_stats.txt"), optional: true, emit: readBinsStats

    shell:
    '''
    echo "test"
    INDEX=!{sample}.index
    bowtie2-build --threads 28 --quiet !{contigs} $INDEX
    bowtie2 -p !{task.cpus} --very-sensitive -x $INDEX --un-conc !{sample}_discordand.fq --interleaved fastq.fq.gz 2> !{sample}_bowtie_stats.txt > out.sam
    paste !{sample}_discordand.1.fq !{sample}_discordand.2.fq | paste - - - - | awk -v OFS="\n" -v FS="\t" '{print($1,$3,$5,$7,$2,$4,$6,$8)}' | pigz --best --processes !{task.cpus} > interleaved.discordand.fq.gz
    '''
}


workflow wUnmappedReadsFile {
   take:
     sampleReadsFile
     sampleBinsFile
   main:
     sampleBinsFile | splitCsv(sep: '\t', header: true) | map { bin -> [bin.SAMPLE, "myMethod",  bin.BINS] } | set{ sampleBinsList } 
     sampleReadsFile | splitCsv(sep: '\t', header: true) \
       | map { sample -> [sample.SAMPLE, sample.READS] } | set {sampleReadsList}     

     _wUnmappedReads(sampleReadsList, sampleBinsList)
   emit:
     unmappedReads = _wUnmappedReads.out.unmappedReads 
     concatSampleBins = _wUnmappedReads.out.concatSampleBins
}


workflow wUnmappedReadsList {
   take:
     sampleReads
     sampleBins 
   main:
     sampleBins | map { sample -> flattenBins(sample) } | flatMap { bin -> bin } | set { flattenedSampleBins }
     _wUnmappedReads(sampleReads, flattenedSampleBins)
   emit:
     unmappedReads = _wUnmappedReads.out.unmappedReads 
     concatSampleBins = _wUnmappedReads.out.concatSampleBins
}


workflow _wUnmappedReads {
   take:
     sampleReads
     sampleBins 
   main:
     BIN_FILE_IDX = 1
     COMBINED_BINS_TMP = params.tempdir + "/combinedBins"
     file(COMBINED_BINS_TMP).mkdirs()
     sampleBins | collectFile(tempDir: COMBINED_BINS_TMP){ sample -> [ "${sample[0]}", file(sample[BIN_FILE_IDX]).text] } \
       | map { bin -> [ bin.baseName, bin ] } | set { concatSampleBins }

     SAMPLE_NAME_IDX = 0
     METHOD_NAME_IDX = 1
     sampleBins | map { sample -> [sample[SAMPLE_NAME_IDX], sample[METHOD_NAME_IDX]] } \
       | unique() | set { sampleMethod }
  
     sampleMethod | join(concatSampleBins, by:SAMPLE_NAME_IDX) | join(sampleReads, by: 0) | pBowtie
   emit:
     unmappedReads = pBowtie.out.unmappedReads
     concatSampleBins = concatSampleBins
}
