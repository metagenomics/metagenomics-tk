nextflow.enable.dsl=2


include { pMinimap2Index as pMinimap2IndexLong; \
          pMapMinimap2 as pMapMinimap2Long; } from './processes'

include { pCovermGenomeCoverage; } from '../binning/processes'

def getModulePath(module){
    return module.name + '/' + module.version.major + "." +
          module.version.minor + "." +
          module.version.patch
}

def getOutput(SAMPLE, RUNID, TOOL, filename){
    return SAMPLE + '/' + RUNID + '/' + params.modules.readMapping.name + '/' + 
         params.modules.readMapping.version.major + "." + 
         params.modules.readMapping.version.minor + "." + 
         params.modules.readMapping.version.patch +
         '/' + TOOL + '/' + filename
}

process pBwaIndex {
    container "${params.bwa_image}"
    label 'large'
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


process pMapBwa {
    label 'large'
    container "${params.samtools_bwa_image}"
    when params.steps.containsKey("readMapping")
    publishDir params.output, mode: "${params.publishDirMode}", saveAs: { filename -> getOutput("${sampleID}", params.runid ,"bwa", filename) }
    input:
      tuple val(sampleID), path(sample, stageAs: "sample*"), path(representatives), path(index, stageAs: "*") 
    output:
      tuple val("${sampleID}"), path("*bam"), path("*bai"), emit: alignment
      tuple file(".command.sh"), file(".command.out"), file(".command.err"), file(".command.log"), emit: logs
    shell:
    template('bwa.sh')
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

     _wReadMappingBwa(ontReads, ontMedianQuality, paired, single, genomesList)
   emit:
     trimmedMean = _wReadMappingBwa.out.trimmedMean
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
     trimmedMean = _wReadMappingBwa.out.trimmedMean
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

     // Combine index with samples
     samplesPaired | join(samplesSingle, remainder: true) \
       | map { sample -> [sample[SAMPLE_NAME_IDX], sample.findAll().tail()] } \
       | combine(shortPairedReadIndex) | set {illumina}

     pMapBwa(illumina)

     // Map ONT data
     samplesONT | combine(ontIndex) | set {ont}
     pMapMinimap2Long(Channel.value(params?.steps.containsKey("readMapping") \
	&& Channel.value(params?.steps?.readMapping.containsKey("minimap"))), ont)
 
     DO_NOT_ESTIMATE_IDENTITY = "-1"
     pMapBwa.out.alignment | combine(genomes | map {it -> file(it)} \
      | toList() | map { it -> [it]})  \
      | combine(Channel.value(DO_NOT_ESTIMATE_IDENTITY)) \
      | set { covermBWAInput }

     pMapMinimap2Long.out.alignment | combine(genomes | map {it -> file(it)} \
      | toList() | map { it -> [it]}) \
      | join(ontMedianQuality, by: SAMPLE_NAME_IDX) | set { covermMinimapInput }

     ALIGNMENT_INDEX = 2
     pCovermGenomeCoverage(Channel.value(params.steps?.readMapping?.find{ it.key == "coverm" }?.value), \
	Channel.value([getModulePath(params.modules.readMapping), \
	"genomeCoverage", params.steps?.readMapping?.coverm?.additionalParams]), \
	covermBWAInput | mix(covermMinimapInput))

   emit:
     trimmedMean = pCovermGenomeCoverage.out.trimmedMean
}
