include { wSaveSettingsList } from '../config/module'


include { pMinimap2Index as pMinimap2IndexLong; \
          pMapMinimap2 as pMapMinimap2Long; } from './processes'

include { pDumpLogs } from '../utils/processes'

include { pCovermGenomeCoverage; } from '../binning/processes'

def getModulePath(module){
    return module.name + '/' + module.version.major + "." +
          module.version.minor + "." +
          module.version.patch
}

def getOutput(RUNID, TOOL, filename){
    return 'AGGREGATED/' + RUNID + '/' + params.modules.readMapping.name + '/' + 
         params.modules.readMapping.version.major + "." + 
         params.modules.readMapping.version.minor + "." + 
         params.modules.readMapping.version.patch +
         '/' + TOOL + '/' + filename
}

process pVerticalConcatFinal {

    label 'tiny'

    cache 'deep'

    publishDir params.output, mode: "${params.publishDirMode}", saveAs: { filename -> getOutput(params.runid, "abundanceMatrix", "${name}") }

    when params.steps.containsKey("readMapping")

    input:
    file('sample?')
    val(name)

    output:
    path("abundance.tsv")

    shell:
    template("verticalConcat.sh")
}


process pVerticalConcat {

    label 'tiny'

    cache 'deep'

    when params.steps.containsKey("readMapping")

    input:
    file('sample?')

    output:
    path("abundance.tsv")

    shell:
    template("verticalConcat.sh")
}





process pBwaIndex {
    container "${params.bwa_image}"
    label 'highmemLarge'
    when params.steps.containsKey("readMapping") && params.steps.readMapping.containsKey("bwa")
    input:
      path(representatives)
    output:
      path('*.{amb,ann,bwt,pac,sa,fa}')
    shell:
      """
      bwa index !{params.steps.readMapping.bwa.additionalParams.bwa_index} !{representatives}
      """
}

process pBwa2Index {
    container "${params.bwa2_image}"
    label 'highmemLarge'
    when params.steps.containsKey("readMapping") && params.steps.readMapping.containsKey("bwa2")
    input:
      path(representatives)
    output:
      path('*.{0123,amb,ann,64,pac}')
    shell:
      """
      bwa-mem2 index !{params.steps.readMapping.bwa2.additionalParams.bwa2_index} !{representatives}
      """
}

process pMapBwa {
    label 'highmemLarge'
    container "${params.samtools_bwa_image}"
    when params.steps.containsKey("readMapping") && params.steps.readMapping.containsKey("bwa")
    publishDir params.output, mode: "${params.publishDirMode}", saveAs: { filename -> getOutput(params.runid ,"bwa", filename) }
    input:
      tuple val(sampleID), path(sample, stageAs: "sample*"), path(representatives), path(index, stageAs: "*") 
    output:
      tuple val("${sampleID}"), path("*bam"), path("*bai"), emit: alignment
      tuple val("${sampleID}"), val("${output}"), val(params.LOG_LEVELS.INFO), file(".command.sh"), \
        file(".command.out"), file(".command.err"), file(".command.log"), emit: logs
    shell:
    output = getOutput(params.runid, "bwa", "")
    template('bwa.sh')
}


process pMapBwa2 {
    label 'highmemLarge'
    container "${params.samtools_bwa2_image}"
    when params.steps.containsKey("readMapping") && params.steps.readMapping.containsKey("bwa2")
    publishDir params.output, mode: "${params.publishDirMode}", saveAs: { filename -> getOutput(params.runid ,"bwa2", filename) }
    input:
      tuple val(sampleID), path(sample, stageAs: "sample*"), path(representatives), path(index, stageAs: "*") 

    output:
      tuple val("${sampleID}"), path("*bam"), path("*bai"), emit: alignment
      tuple val("${sampleID}"), val("${output}"), val(params.LOG_LEVELS.INFO), file(".command.sh"), \
        file(".command.out"), file(".command.err"), file(".command.log"), emit: logs

    shell:
    output = getOutput(params.runid, "bwa2", "")
    template('bwa2.sh')
}


/*
*
* wListReadMappingBwa maps a set of fastq samples against a set of genomes.
* Input Parameters
*
*  * ONT, paired and single end sample sheets must have the columns "READS" and "SAMPLE"  
*  * MAGs file contains paths to every input genome (BINs column).
*  * Columns of the mags file must be SAMPLE for sample identifier and READS for the path to the fastq files.
*
*/
workflow wFileReadMappingBwa {
   main:
    
     paired = Channel.empty()
     single = Channel.empty()
     ontReads = Channel.empty()
     ontMedianQuality = Channel.empty()
   
     GENOMES_PATH_IDX = 0
     Channel.from(file(params?.steps?.readMapping?.mags)) \
       | splitCsv(sep: '\t', header: false, skip: 1) \
       | map { it -> file(it[GENOMES_PATH_IDX]) } \
       | set {genomesList}

     if(params?.steps?.readMapping?.samples.containsKey("paired")){
       Channel.from(file(params?.steps?.readMapping?.samples?.paired)) \
         | splitCsv(sep: '\t', header: true)\
         | map { it -> [it.SAMPLE, it.READS] } | set {paired}
     }     

     if(params?.steps?.readMapping?.samples.containsKey("single")){
       Channel.from(file(params?.steps?.readMapping?.samples?.single)) \
         | splitCsv(sep: '\t', header: true)\
         | map { it -> [it.SAMPLE, it.READS] } | set {single}
     }

     if(params?.steps?.readMapping?.samples.containsKey("ont")){
       Channel.from(file(params?.steps?.readMapping?.samples?.ont)) \
         | splitCsv(sep: '\t', header: true) | set { ont }
       ont | map { sample -> [sample.SAMPLE, sample.READS] } | set {ontReads}
       ont | map { sample -> [sample.SAMPLE, sample.MEDIAN_PHRED] } | set {ontMedianQuality}
     }   

     SAMPLE_IDX = 0
     wSaveSettingsList(paired | mix(single) \
	| mix(ontReads) | map { it -> it[SAMPLE_IDX] })

     _wReadMappingBwa(ontReads, ontMedianQuality, paired, single, genomesList)
}

/*
*
* wListReadMappingBwa maps a set of fastq samples against a set of genomes.
* Input Parameters
*
*  * Entries in the samples channels (ont, paired and single illumina) consist of lists with the format [SAMPLE, READS]. 
*
*  * Entries in the ontMedianQuality channel must have the format [SAMPLE, QUALITY].
*   
*  * Entries of the genomes channel are paths to the genomes. 
*
*/
workflow wListReadMappingBwa {
   take:
     samplesONT
     ontMedianQuality
     samplesPaired
     samplesSingle
     genomes
   main:
     _wReadMappingBwa(samplesONT, ontMedianQuality, samplesPaired, samplesSingle, genomes)
   emit:
     trimmedMeanMatrix = _wReadMappingBwa.out.trimmedMeanMatrix
     relativeAbundanceMatrix = _wReadMappingBwa.out.relativeAbundanceMatrix
}

workflow _wCreateTrimmedMeanMatrix {
  take:
    countColumn
  main:
    _wCreateMatrix(countColumn, Channel.value("trimmedMeanMatrix.tsv"))
  emit:
    trimmedMeanMatrix = _wCreateMatrix.out.abundanceMatrix
}

workflow _wCreateRelativeAbundanceMatrix {
  take:
    countColumn
  main:
    _wCreateMatrix(countColumn, Channel.value("relativeAbundanceMatrix.tsv"))
  emit:
    relativeAbundanceMatrix = _wCreateMatrix.out.abundanceMatrix
}

workflow _wCreateRPKMMatrix {
  take:
    countColumn
  main:
    _wCreateMatrix(countColumn, Channel.value("rpkmMatrix.tsv"))
  emit:
    relativeAbundanceMatrix = _wCreateMatrix.out.abundanceMatrix
}

workflow _wCreateTPMMatrix {
  take:
    countColumn
  main:
    _wCreateMatrix(countColumn, Channel.value("tpmMatrix.tsv"))
  emit:
    relativeAbundanceMatrix = _wCreateMatrix.out.abundanceMatrix
}


workflow _wCreateMatrix {
   take:
     countColumn
     name
   main:
     COLLECT_BUFFER=10000

     // Collect abundance/ count information of every sample
     countColumn | combine(countColumn | count )  \
	| filter({ sample, abundance, count -> count > 1 }) \
	| map { sample, abundance, count -> file(abundance) }\
        | buffer( size: COLLECT_BUFFER, remainder: true ) \
        | pVerticalConcat | collect | set { chunks }
     pVerticalConcatFinal(chunks, name) | set { abundanceMatrix }
    emit:
      abundanceMatrix
}


workflow _wReadMappingBwa {
   take:
     samplesONT
     ontMedianQuality
     samplesPaired
     samplesSingle
     genomes
   main:
     SAMPLE_NAME_IDX=0

     // Create temporary directory for merged fasta files
     genomesTempDir = params.tempdir + "/genomesToBeMerged"
     file(genomesTempDir).mkdirs()

     // Concatenate genomes 
     genomes | collectFile(tempDir: genomesTempDir){ item -> [ "mergedGenomes.fasta", item.text ] } \
      | set { genomesMerged }

     // Create BWA and Minimap index of all genomes
     BWA_INDEX_IDX=0
     GENOMES_IDX=1
     MERGED_GENOME_IDX=2
      
     pMinimap2IndexLong(Channel.value(params?.steps.containsKey("readMapping") \
	&& Channel.value(params?.steps?.readMapping.containsKey("minimap"))), \
	Channel.value("map-ont"), samplesONT | combine(genomesMerged) \
	| map { sample -> sample[MERGED_GENOME_IDX] } ) | set { ontIndex }

     // Create BWA index
     genomesMerged | pBwaIndex | map{ bwaIndex -> [bwaIndex]} \
      | combine(genomesMerged)  \
      | map{ it -> [it[GENOMES_IDX], it[BWA_INDEX_IDX]] } \
      | set  {shortPairedReadIndex}

     // Create BWA2 index
     genomesMerged | pBwa2Index | map{ bwa2Index -> [bwa2Index]} \
      | combine(genomesMerged)  \
      | map{ it -> [it[GENOMES_IDX], it[BWA_INDEX_IDX]] } \
      | set  {shortPairedReadIndex}

     // Combine index with samples
     samplesPaired | join(samplesSingle, remainder: true) \
       | map { sample -> [sample[SAMPLE_NAME_IDX], sample.findAll().tail()] } \
       | combine(shortPairedReadIndex) | set {illumina}

     pMapBwa(illumina)

     pMapBwa2(illumina)

     // Map ONT data
     samplesONT | combine(ontIndex) | set {ont}
     pMapMinimap2Long(Channel.value(params?.steps.containsKey("readMapping") \
	&& Channel.value(params?.steps?.readMapping.containsKey("minimap"))), ont)
 
     DO_NOT_ESTIMATE_IDENTITY = "-1"
     pMapBwa.out.alignment | mix(pMapBwa2.out.alignment ) | combine(genomes | map {it -> file(it)} \
      | toList() | map { it -> [it]})  \
      | combine(Channel.value(DO_NOT_ESTIMATE_IDENTITY)) \
      | set { covermBWAInput }

     pMapMinimap2Long.out.alignment | combine(genomes | map {it -> file(it)} \
      | toList() | map { it -> [it]}) \
      | join(ontMedianQuality, by: SAMPLE_NAME_IDX) | set { covermMinimapInput }

     ALIGNMENT_INDEX = 2
     pCovermGenomeCoverage(Channel.value(params.steps?.readMapping?.find{ it.key == "coverm" }?.value), \
	Channel.value("AGGREGATED"), \
	Channel.value([getModulePath(params.modules.readMapping), \
	"genomeCoverage", params.steps?.readMapping?.coverm?.additionalParams]), \
	covermBWAInput | mix(covermMinimapInput))

     pMapBwa2.out.logs | mix(pMapBwa.out.logs) \
	| mix(pMapMinimap2Long.out.logs) | mix(pCovermGenomeCoverage.out.logs) | pDumpLogs

     _wCreateTrimmedMeanMatrix(pCovermGenomeCoverage.out.trimmedMean)
     _wCreateRelativeAbundanceMatrix(pCovermGenomeCoverage.out.relativeAbundance)
     _wCreateRPKMMatrix(pCovermGenomeCoverage.out.rpkm)
     _wCreateTPMMatrix(pCovermGenomeCoverage.out.tpm)
   emit:
     trimmedMeanMatrix = _wCreateTrimmedMeanMatrix.out.trimmedMeanMatrix
     relativeAbundanceMatrix = _wCreateRelativeAbundanceMatrix.out.relativeAbundanceMatrix
}
