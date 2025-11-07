include { wSaveSettingsList } from '../config/module'

include { pBowtie2; pMinimap2; pBwa; pBwa2; pGetBinStatistics; \
	pCovermGenomeCoverage; pCovermContigsCoverage; } from  '../binning/processes'

include { pCovermCount as pCovermCount; pCovermCount as pCovermCountONT; } from './processes'

include { pMashPaste as pMashPasteChunk; \
	  pMashPaste as pMashPasteFinal; } from  '../dereplication/bottomUpClustering/processes'
include { pDumpLogs } from '../utils/processes'



def getOutput(SAMPLE, RUNID, TOOL, filename){
    return SAMPLE + '/' + RUNID + '/' + params.modules.fragmentRecruitment.name + '/' + 
         params.modules.fragmentRecruitment.version.major + "." +
         params.modules.fragmentRecruitment.version.minor + "." +
         params.modules.fragmentRecruitment.version.patch +
         '/' + TOOL + '/' + filename
}


process pMashScreen {

    label 'highmemMedium'

    tag "Sample: $sample"

    cache 'deep'

    publishDir params.output, mode: "${params.publishDirMode}", saveAs: { filename -> getOutput("${sample}", params.runid, "mashScreen",  filename) }

    container "${params.mash_image}"

    containerOptions (params.apptainer ? "" : Utils.getDockerNetwork())

    when params?.steps?.fragmentRecruitment?.mashScreen != null

    input:
    tuple val(sample), file(pairedReads), file(singleReads), file(sketch) 

    output:
    tuple val("${sample}"), file("mash_screen.tsv"), optional: true, emit: mashScreenOutput
    tuple val("${sample}"), file("selected_genomes.tsv"), optional: true, emit: mashScreenFilteredOutput
    tuple file(".command.sh"), file(".command.out"), file(".command.err"), file(".command.log")

    shell:
    template "mashScreen.sh"
}


process pGenomeContigMapping {

    container "${params.ubuntu_image}"

    tag "Sample: $sample"

    label 'tiny'

    publishDir params.output, mode: "${params.publishDirMode}", saveAs: { filename -> getOutput("${sample}", params.runid, "stats", filename) }

    input:
    tuple val(sample), path(genomes)

    output:
    tuple val("${sample}"), path("${sample}_genome_contig_mapping.tsv"), emit: mapping
    tuple file(".command.sh"), file(".command.out"), file(".command.err"), file(".command.log"), emit: logs

    shell:
    '''
    for g in $(ls -1 !{genomes}); do
       seqkit fx2tab -i -n ${g} | sed "s/^/${g}\t/g" >> !{sample}_genome_contig_mapping.tsv
    done
    '''
}



process pSaveMatchedGenomes {

    container "${params.ubuntu_image}"

    tag "Sample: $sample"

    label 'tiny'

    publishDir params.output, mode: "${params.publishDirMode}", \
	saveAs: { filename -> getOutput("${sample}", params.runid, "", filename) }, \
        pattern: "{matches/*}"

    input:
    tuple val(sample), path(genomes)

    output:
    path("matches/*"), emit: matches
    tuple val("${sample}"), val("${output}"), val(params.LOG_LEVELS.ALL), file(".command.sh"), \
      file(".command.out"), file(".command.err"), file(".command.log"), emit: logs

    shell:
    output = getOutput("${sample}", params.runid, "matches", "")
    '''
    mkdir matches
    cp !{genomes} matches
    '''
}



/*
*
* This workflow takes paired end reads as input and searches for 
* possible matches in all user provided genomes.
* The input reads tsv file must have the columns: SAMPLE, READS
* The input genomes tsv file must have the columns: PATH 
*
*/
workflow wMashScreenFile {
   main:
     pairedReads = Channel.empty()
     if(params.steps.fragmentRecruitment.containsKey("mashScreen")){
        Channel.from(file(params.steps.fragmentRecruitment.mashScreen.samples.paired)) | splitCsv(sep: '\t', header: true) \
             | map { line -> [ line.SAMPLE, file(line.READS)]} | set { pairedReads  }
     }

     singleReads = Channel.empty()
     if(params.steps.fragmentRecruitment.containsKey("mashScreen") && params.steps.fragmentRecruitment.mashScreen.samples.containsKey("single")){
     	Channel.from(file(params.steps.fragmentRecruitment.mashScreen.samples.single)) | splitCsv(sep: '\t', header: true) \
             | map { line -> [ line.SAMPLE, file(line.READS)]} | set { singleReads  }
     }

     ontReads = Channel.empty()
     if(params.steps.fragmentRecruitment.containsKey("mashScreen") && params.steps.fragmentRecruitment.mashScreen.samples.containsKey("ont")){
     	Channel.from(file(params.steps.fragmentRecruitment.mashScreen.samples.ont)) | splitCsv(sep: '\t', header: true) \
             | map { line -> [ line.SAMPLE, file(line.READS)]} | set { ontReads  }
     }
 
     SAMPLE_IDX = 0
     wSaveSettingsList(ontReads | mix(singleReads) | mix(pairedReads) | map { it -> it[SAMPLE_IDX] } | unique)

     _wMashScreen(pairedReads, singleReads, ontReads, Channel.empty())
}


/*
*
* This workflow takes paired end reads as input and searches for 
* possible matches in all user provided genomes.
* The input must have the form: [SAMPLEID, Interleaved reads]
*
*/
workflow wMashScreenList {
   take:
      pairedReads   
      ontReads
      medianQuality
   main:
     _wMashScreen(pairedReads, Channel.empty(), ontReads, medianQuality)

   emit:
     foundGenomesPerSample = _wMashScreen.out.foundGenomesPerSample
     foundGenomesSeperated = _wMashScreen.out.foundGenomesSeperated
     binsStats = _wMashScreen.out.binsStats
     contigCoverage = _wMashScreen.out.contigCoverage
     genomeContigMapping = _wMashScreen.out.genomeContigMapping
}


/*
*
* Method takes two channels with map entries and two keys as input.
* Channels are joined by the keys provided.
* Resulting channel is returned as output.
*
*/
def mapJoin(channel_a, channel_b, key_a, key_b){
    channel_a \
        | map{ it -> [it[key_a], it] } \
        | cross(channel_b | map{it -> [it[key_b], it]}) \
        | map { it[0][1] + it[1][1] }
}


workflow _wRunMash {
   take:
     sampleReads
     singleReads
     genomes
   main:
     BUFFER_SKETCH = 1000

     pMashSketchGenomeGroup(params?.steps.containsKey("fragmentRecruitment") &&  params?.steps.fragmentRecruitment.containsKey("mashScreen"), \
	Channel.value(params.steps?.fragmentRecruitment?.mashScreen?.additionalParams?.mashSketch), genomes | buffer(size: BUFFER_SKETCH, remainder: true))

     // Combine sketches
     pMashSketchGenomeGroup.out.sketches | flatten | toSortedList | flatten | buffer(size: BUFFER_SKETCH, remainder: true) | set { mashPasteInput }

     pMashPasteChunk(params?.steps.containsKey("fragmentRecruitment") &&  params?.steps.fragmentRecruitment.containsKey("mashScreen"), \
	Channel.value([Utils.getModulePath(params.modules.fragmentRecruitment), "mash/paste"]),  mashPasteInput)

     pMashPasteFinal(params?.steps.containsKey("fragmentRecruitment") &&  params?.steps.fragmentRecruitment.containsKey("mashScreen"), \
	Channel.value([Utils.getModulePath(params.modules.fragmentRecruitment), "mash/paste"]),  pMashPasteChunk.out.sketch | collect(flat: false))

     pMashPasteFinal.out.logs | pDumpLogs

     // Screen reads for genomes
     SAMPLE_IDX = 0
     sampleReads | join(singleReads, by: SAMPLE_IDX, remainder: true) | combine(pMashPasteFinal.out.sketch) \
    	| pMashScreen
   emit:
     mashScreenFilteredOutput = pMashScreen.out.mashScreenFilteredOutput 
}


workflow _wRunMapping {
  take:
     mashOutput
     sampleReads
     singleReads
     ontReads
     genomesMap
  main:
     SAMPLE_IDX = 0

     fragmentRecruitmentGenomes = params.tempdir + "/fragmentRecruitmentGenomes"
     file(fragmentRecruitmentGenomes).mkdirs()

     // Get found genomes and merge them as a preparation for short read mapper 
     GENOMES_IDX = 1
     FIRST_GENOME_IDX = 0
     SAMPLE_2_IDX = 1
     PATH_2_IDX = 2
     mashOutput | splitCsv(sep: '\t', header: false) \
	| map { line -> [line[GENOMES_IDX][FIRST_GENOME_IDX], line[SAMPLE_IDX]] } \
	| combine(genomesMap, by: SAMPLE_IDX) \
	| collectFile(tempDir: fragmentRecruitmentGenomes){ item -> [ "${item[SAMPLE_2_IDX]}", item[PATH_2_IDX].text ] } \
	| map { genome -> [genome.name, genome] } | set{genomesMerged}

     // Validate if the found genomes via mash can also be detected via a read mapper and coverm
     emptyFile = file(params.tempdir + "/empty")
     emptyFile.text = ""
     GENOMES_MERGED_IDX = 2
     GENOMES_MERGED_2_IDX = 1
     SINGLE_SAMPLE_IDX = 3
     SAMPLE_3_IDX = 2
     sampleReads | join(genomesMerged, by: SAMPLE_IDX) \
	| map { sample -> [sample[SAMPLE_IDX], sample[GENOMES_MERGED_IDX], sample[SAMPLE_2_IDX]] } \
	| join(singleReads, by: SAMPLE_IDX, remainder: true) \
	| filter { sample -> sample[SAMPLE_2_IDX] != null } \
	| map { sample -> sample[SINGLE_SAMPLE_IDX]? sample : [sample[SAMPLE_IDX], sample[GENOMES_MERGED_2_IDX], sample[SAMPLE_3_IDX], emptyFile] } \
	| set { mapperInput }

     pBowtie2(Channel.value(params?.steps.containsKey("fragmentRecruitment")  \
 	&& params.steps?.fragmentRecruitment?.mashScreen?.additionalParams.containsKey("bowtie")), \
	Channel.value([Utils.getModulePath(params.modules.fragmentRecruitment), \
        "readMapping/bowtie", params.steps?.fragmentRecruitment?.mashScreen?.additionalParams?.bowtie, \
	params.steps?.fragmentRecruitment?.mashScreen?.additionalParams?.samtoolsViewBowtie, \
	params.steps.containsKey("fragmentRecruitment")]), mapperInput)

     pBwa(Channel.value(params?.steps.containsKey("fragmentRecruitment") \
	&& params.steps?.fragmentRecruitment?.mashScreen?.additionalParams.containsKey("bwa")), \
	Channel.value([Utils.getModulePath(params.modules.fragmentRecruitment), \
        "readMapping/bwa", params.steps?.fragmentRecruitment?.mashScreen?.additionalParams?.bwa, \
	params.steps?.fragmentRecruitment?.mashScreen?.additionalParams?.samtoolsViewBwa, \
	params.steps.containsKey("fragmentRecruitment")]), mapperInput)

     pBwa2(Channel.value(params?.steps.containsKey("fragmentRecruitment") \
	&& params.steps?.fragmentRecruitment?.mashScreen?.additionalParams.containsKey("bwa2")), \
	Channel.value([Utils.getModulePath(params.modules.fragmentRecruitment), \
        "readMapping/bwa2", params.steps?.fragmentRecruitment?.mashScreen?.additionalParams?.bwa2, \
	params.steps?.fragmentRecruitment?.mashScreen?.additionalParams?.samtoolsViewBwa2, \
	params.steps.containsKey("fragmentRecruitment")]), mapperInput)

     ontReads | join(genomesMerged, by: SAMPLE_IDX) \
	| map { sample -> [sample[SAMPLE_IDX], sample[GENOMES_MERGED_IDX], sample[SAMPLE_2_IDX]] } \
	| set { minimapInput }

     pMinimap2(Channel.value(params?.steps.containsKey("fragmentRecruitment")), \
	Channel.value([Utils.getModulePath(params.modules.fragmentRecruitment), \
        "readMapping/minimap",  params.steps?.fragmentRecruitment?.mashScreen?.additionalParams?.minimap, \
        params.steps?.fragmentRecruitment?.mashScreen?.additionalParams?.samtoolsViewMinimap, \
	params.steps.containsKey("fragmentRecruitment")]), minimapInput)

      pBowtie2.out.mappedReads | mix(pBwa.out.mappedReads, pBwa2.out.mappedReads) | set { mappedReads }

  emit:
    mappedShortReads = mappedReads
    minimapMappedReads = pMinimap2.out.mappedReads
       
}


workflow _wGetStatistics {
  take:
    mashScreenFilteredOutput
    mappedShortReads
    minimapMappedReads
    ontMedianQuality
    genomesMap
  main:
     SAMPLE_IDX = 0

     GENOMES_IDX = 1
     GENOMES_2_IDX = 2
     FIRST_GENOME_IDX = 0
     SAMPLE_2_IDX = 1
     PATH_2_IDX = 2

     DO_NOT_SET_IDENTITY_AUTOMATICALLY = "-1"

     // Create input for coverm
     mashScreenFilteredOutput | splitCsv(sep: '\t', header: false) \
	| map { line -> [line[GENOMES_IDX][FIRST_GENOME_IDX], line[SAMPLE_IDX]] } \
	| combine(genomesMap, by: SAMPLE_IDX) | map { sample -> [sample[SAMPLE_2_IDX], sample[PATH_2_IDX]] }  \
	| groupTuple(by: SAMPLE_IDX) | set { covermGenomesInput }

     mappedShortReads \
	| join(covermGenomesInput, by: SAMPLE_IDX) \
	| combine(Channel.value(DO_NOT_SET_IDENTITY_AUTOMATICALLY)) \
        | set { covermBowtieReadsInput  }

     minimapMappedReads \
	| join(covermGenomesInput, by: SAMPLE_IDX) \
	| join(ontMedianQuality, by: SAMPLE_IDX) \
        | set { covermMinimapReadsInput }

     pCovermCount(Channel.value(params?.steps?.fragmentRecruitment?.mashScreen?.additionalParams.find{ it.key == "coverm"}?.value), \
	covermBowtieReadsInput)

     pCovermCountONT(Channel.value(params?.steps?.fragmentRecruitment?.mashScreen?.additionalParams.find{ it.key == "covermONT"}?.value), \
	covermMinimapReadsInput)

     pCovermCount.out.foundGenomes | mix(pCovermCountONT.out.foundGenomes) | set { foundGenomes } 

     // Found genomes reported by coverm are just file names. These lines join unique file names with file paths. 
     foundGenomes |  splitCsv(sep: '\t', header: false) \
	| map { genomes -> genomes[GENOMES_IDX] } | flatten | join(genomesMap) | map { genome -> genome[GENOMES_IDX] } \
	| set { covermFilteredGenomes }

     // Get Bin coverage statistics of the alignment
     covermFilteredGenomes | map { genome -> ["EXTERNAL_GENOMES", file(genome).name, genome] } | set { genomesSeperated }

     foundGenomes | splitCsv(sep: '\t', header: false) \
        | map { genomes -> Utils.flattenTuple(genomes).flatten()} \
	| combine(genomesMap | map { genome -> genome.reverse() }, by: GENOMES_IDX) \
	| map { genomes -> [genomes[SAMPLE_2_IDX], genomes[GENOMES_2_IDX]]} \
	| set { foundGenomesSeperated }

     foundGenomesSeperated | groupTuple(by: SAMPLE_IDX) | set { foundGenomesInGroup }
     foundGenomesSeperated | map { genome -> [genome[SAMPLE_IDX], file(genome[GENOMES_IDX]).name, genome[GENOMES_IDX]] } \
	| set { foundGenomesIDSeperated }
  
     foundGenomesInGroup | pGenomeContigMapping

     foundGenomesInGroup | pSaveMatchedGenomes


     pGenomeContigMapping.out.mapping | join(mappedShortReads, by: SAMPLE_IDX) \
        | combine(Channel.from("stats")) | join(foundGenomesInGroup, by: SAMPLE_IDX) \
        | combine(Channel.value(DO_NOT_SET_IDENTITY_AUTOMATICALLY)) \
        | set { shortReadMappingStatsInput }

     pGenomeContigMapping.out.mapping | join(minimapMappedReads, by: SAMPLE_IDX) \
        | combine(Channel.from("external")) | join(foundGenomesInGroup, by: SAMPLE_IDX) \
        | join(ontMedianQuality, by: SAMPLE_IDX)
        | set { minimapMappingStatsInput }

     pGetBinStatistics(Utils.getModulePath(params.modules.fragmentRecruitment), shortReadMappingStatsInput | mix(minimapMappingStatsInput))

     pCovermContigsCoverage(Channel.value(params?.steps?.fragmentRecruitment.find{ it.key == "contigsCoverage"}?.value), \
	Channel.value([Utils.getModulePath(params.modules.fragmentRecruitment), \
	"contigCoverage", params?.steps?.fragmentRecruitment?.contigsCoverage?.additionalParams]), \
	mappedShortReads | combine(Channel.value(DO_NOT_SET_IDENTITY_AUTOMATICALLY)) \
	| mix(minimapMappedReads | join(ontMedianQuality, by: SAMPLE_IDX)))

     emptyFile = file(params.tempdir + "/empty")
     emptyFile.text = ""

     minimapMappedReads | join(foundGenomesInGroup, by: SAMPLE_IDX) \
        | join(ontMedianQuality, by: SAMPLE_IDX) | set { minimapMappedReadsCovInput }

     mappedShortReads | join(foundGenomesInGroup, by: SAMPLE_IDX) \
        | combine(Channel.value(DO_NOT_SET_IDENTITY_AUTOMATICALLY)) | set { mappedShortReadsCovInput }

     ALIGNMENT_INDEX = 2
     pCovermGenomeCoverage(Channel.value(params?.steps?.fragmentRecruitment.find{ it.key == "genomeCoverage"}?.value), \
        Channel.value(""), \
	Channel.value([Utils.getModulePath(params.modules.fragmentRecruitment), \
	"genomeCoverage", params?.steps?.fragmentRecruitment?.genomeCoverage?.additionalParams]), \
	 minimapMappedReadsCovInput | mix(mappedShortReadsCovInput) \
	| map { sample -> sample.addAll(ALIGNMENT_INDEX, emptyFile); sample })

     pSaveMatchedGenomes.out.logs | mix(pCovermGenomeCoverage.out.logs) | pDumpLogs

  emit:
     foundGenomesPerSample = foundGenomesInGroup 
     foundGenomesSeperated = foundGenomesIDSeperated
     genomesSeperated = genomesSeperated
     genomeContigMapping = pGenomeContigMapping.out.mapping 
     binsStats = pGetBinStatistics.out.binsStats     
     contigCoverage = pCovermContigsCoverage.out.coverage
}


process pUnzipGroup {

  label 'tiny'

  fair true

  container "${params.ubuntu_image}"

  when params?.steps.containsKey("fragmentRecruitment") && params.steps.fragmentRecruitment.containsKey("mashScreen")

  time params?.steps.containsKey("fragmentRecruitment") \
	&& params.steps.fragmentRecruitment.containsKey("mashScreen") \
	? Utils.setTimeLimit(params.steps.fragmentRecruitment.mashScreen.unzip, params.modules.fragmentRecruitment.process.unzip.defaults, params.resources.tiny) : ""

  cache 'deep'

  input:
  path(x, stageAs: "input/*")

  output:
  path("*", type: "file")

  shell:
  '''
  for f in input/*; do  
	cp $f . ; 
	name=$(basename $f); 
	if gzip -t $name; then 
		gunzip  -d $name ;
	fi 
  done
  '''
}


process pMashSketchGenomeGroup {

    container "${params.mash_image}"

    label 'tiny'

    fair true

    when:
    run

    input:
    val(run)
    val(mashSketchParams)
    path(x)

    output:
    path("*.msh"), emit: sketches

    shell:
    '''
    for f in * ; do  
    	mash sketch !{mashSketchParams} $f -o $(basename $f).msh
    done
    '''
}


workflow _wMashScreen {
   take: 
     sampleReads
     singleReads
     ontReads
     ontMedianQuality
   main:
     SAMPLE_IDX = 0
     BIN_ID_IDX = 1
     PATH_2_IDX = 2
     STATS_PATH = 1

     genomes = Channel.empty()
     if(params?.steps.containsKey("fragmentRecruitment") \
	&& params?.steps.fragmentRecruitment.containsKey("mashScreen") \
	&& params?.steps.fragmentRecruitment.mashScreen.containsKey("genomes")){

       BUFFER = 500
       // Some tools can not handle gzipped files. Unzip genomes before fragment recruitment
       Channel.fromPath(params?.steps.fragmentRecruitment?.mashScreen?.genomes) | splitCsv(sep: '\t', header: true) \
         | map { line -> file(line.PATH)} | buffer( size: BUFFER, remainder: true ) | pUnzipGroup \
	 |  flatten | set { genomes }
     }

     UNIQUE_IDX=0
     ALL_IDX=0
     genomes | map { genome -> genome.name } | unique | count | combine(genomes | count) \
	| filter { count -> count[UNIQUE_IDX]!=count[ALL_IDX] } | view { error "Genome file names must be unique!" } 

     genomes | map { genome -> [file(genome).name, genome] } | set { genomesMap }

     // Run mash and return matched genomes per sample
     _wRunMash(sampleReads | mix(ontReads), singleReads, genomes)

     // Concatenate genomes per sample and map reads per sample against concatenated genomes via Bowtie 
     _wRunMapping(_wRunMash.out.mashScreenFilteredOutput, sampleReads, singleReads, ontReads, genomesMap)

     // Calculate different mapping statistics per genome
     _wGetStatistics(_wRunMash.out.mashScreenFilteredOutput, \
	 _wRunMapping.out.mappedShortReads, \
         _wRunMapping.out.minimapMappedReads, \
	ontMedianQuality, genomesMap)

     _wGetStatistics.out.genomesSeperated \
	| map { genome -> [SAMPLE: genome[SAMPLE_IDX], BIN_ID: genome[BIN_ID_IDX], PATH: genome[PATH_2_IDX]] } \
	| set {binMap}

     _wGetStatistics.out.binsStats | map { stats -> file(stats[STATS_PATH]) } \
        | splitCsv(sep: '\t', header: true) | set { binsStats }

     mapJoin(binsStats, binMap, "BIN_ID", "BIN_ID") | set {binMap}

     binMap | unique { bin -> bin.BIN_ID } | set { binMap }
    emit:
     foundGenomesPerSample = _wGetStatistics.out.foundGenomesPerSample
     binsStats = binMap
     genomeContigMapping = _wGetStatistics.out.genomeContigMapping
     contigCoverage = _wGetStatistics.out.contigCoverage
     foundGenomesSeperated = _wGetStatistics.out.foundGenomesSeperated
}
