nextflow.enable.dsl=2

include { pBowtie2; pMinimap2; pBwa; pGetBinStatistics; } from  '../binning/processes'
include { pMashSketchGenome; \
	  pMashPaste as pMashPasteChunk; \
	  pMashPaste as pMashPasteFinal; } from  '../dereplication/pasolli/processes'



def getOutput(SAMPLE, RUNID, TOOL, filename){
    return SAMPLE + '/' + RUNID + '/' + params.modules.fragmentRecruitment.name + '/' + 
         params.modules.fragmentRecruitment.version.major + "." +
         params.modules.fragmentRecruitment.version.minor + "." +
         params.modules.fragmentRecruitment.version.patch +
         '/' + TOOL + '/' + filename
}

def getModulePath(module){
    return module.name + '/' + module.version.major + "." +
          module.version.minor + "." +
          module.version.patch
}


process pMashScreen {

    label 'medium'

    tag "Sample: $sample"

    publishDir params.output, mode: "${params.publishDirMode}", saveAs: { filename -> getOutput("${sample}", params.runid, "mashScreen",  filename) }

    container "${params.mash_image}"

    when params?.steps?.fragmentRecruitment?.mashScreen != null

    input:
    tuple val(sample), file(pairedReads), file(singleReads), file(sketch) 

    output:
    tuple val("${sample}"), file("mash_screen.tsv"), emit: mashScreenOutput
    tuple val("${sample}"), file("selected_genomes.tsv"), emit: mashScreenFilteredOutput
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
     genomes = _wMashScreen.out.genomes
     genomesSeperated = _wMashScreen.out.genomesSeperated
     binsStats = _wMashScreen.out.binsStats
}


process pCovermCount {

    label 'small'

    container "${params.ubuntu_image}"

    publishDir params.output, mode: "${params.publishDirMode}", saveAs: { filename -> getOutput("${sample}", params.runid, "coverm", filename) }

    input:
      tuple val(sample), file(mapping), file(listOfRepresentatives), val(medianQuality)

    output:
      tuple val("${sample}"), path("${sample}_stats_out/coveredBases.tsv"), emit: mean
      tuple val("${sample}"), path("${sample}_stats_out/metrics.tsv"), emit: metrics
      tuple val("${sample}"), path("foundGenomes.tsv"), optional: true, emit: foundGenomes
      tuple file(".command.sh"), file(".command.out"), file(".command.err"), file(".command.log")

    shell:
    DO_NOT_ESTIMATE_QUALITY = -1 
    MEDIAN_QUALITY=Double.parseDouble(medianQuality)
    percentIdentity = MEDIAN_QUALITY != DO_NOT_ESTIMATE_QUALITY ? \
	" --min-read-percent-identity "+Utils.getMappingIdentityParam(MEDIAN_QUALITY) : " "
    '''
    OUT=!{sample}_stats_out
    mkdir $OUT
    readlink -f !{listOfRepresentatives} > list.txt 
   
    # Create a mapping between file path basename of the file without ending. (/path/test.1.tsv --> test.1) 
    paste -d '\t' list.txt  <(cat list.txt  | rev | cut -d '/' -f 1  | cut -d '.' -f 2- | rev) > mapping.tsv
    
    # Get covered bases
    coverm genome -t !{task.cpus} -b !{mapping} \
         !{params.steps?.fragmentRecruitment?.mashScreen?.additionalParams?.coverm}  !{percentIdentity}  \
        --genome-fasta-list list.txt --methods covered_bases --output-file covTmpContent.tsv \

    # Get length
    coverm genome -t !{task.cpus} -b !{mapping} --min-covered-fraction 0  \
        --genome-fasta-list list.txt --methods length --output-file lengthTmpContent.tsv \

    # Join length and covered bases
    join -t$'\t' -1 1 -2 1 covTmpContent.tsv lengthTmpContent.tsv > covLengthTmpContent.tsv

    # Exchange header ad add covered fraction column
    sed -i  -e '1 s/^.*$/SAMPLE\tGENOME\tCOVERED_BASES\tLENGTH/' -e "2,$ s/^/!{sample}\t/g" covLengthTmpContent.tsv  \
                && echo "COVERED_FRACTION" > covLengthTmp.tsv \
		&& awk '(NR>1){ tmp=($3/($4/100)) ; printf"%0.2f\\n", tmp }' covLengthTmpContent.tsv >> covLengthTmp.tsv \
                && paste -d$'\t' covLengthTmpContent.tsv covLengthTmp.tsv > $OUT/coveredBases.tsv || true

    # Run other metrics like RPKM, TPM, ...
    coverm genome  -t !{task.cpus} -b !{mapping} \
         !{params.steps?.fragmentRecruitment?.mashScreen?.additionalParams?.coverm} !{percentIdentity} \
	--genome-fasta-list list.txt --methods mean trimmed_mean variance length count reads_per_base rpkm tpm \
	| sed -e '1 s/^.*$/SAMPLE\tGENOME\tMEAN\tTRIMMED_MEAN\tVARIANCE\tLENGTH\tREAD_COUNT\tREADS_PER_BASE\tRPKM\tTPM/' \
	| sed -e "2,$ s/^/!{sample}\t/g" > $OUT/metrics.tsv || true

    coveredBasesCutoff=!{params.steps?.fragmentRecruitment?.mashScreen?.coveredBasesCutoff}
    FOUND_GENOMES=foundGenomes.tsv
    for b in $(awk -v coveredBases=${coveredBasesCutoff} '(NR>1){if ($5 > coveredBases) print $2}' $OUT/coveredBases.tsv); do 
        file=$(grep -P "\t$b$" mapping.tsv | cut -f 1);
	echo $(basename $file)  >> ${FOUND_GENOMES} 
    done
    '''

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
     genomesMap
   main:
     BUFFER_SKETCH = 1000
     PATH_IDX = 0

     pMashSketchGenome(params?.steps.containsKey("fragmentRecruitment") &&  params?.steps.fragmentRecruitment.containsKey("mashScreen"), \
	Channel.value(params.steps?.fragmentRecruitment?.mashScreen?.additionalParams?.mashSketch), genomesMap)

     // Combine sketches
     pMashSketchGenome.out.sketch  | toSortedList | flatten | buffer(size: BUFFER_SKETCH, remainder: true) | set { mashPasteInput }

     pMashPasteChunk(params?.steps.containsKey("fragmentRecruitment") &&  params?.steps.fragmentRecruitment.containsKey("mashScreen"), \
	Channel.value([getModulePath(params.modules.fragmentRecruitment), "mash/paste"]),  mashPasteInput)

     pMashPasteFinal(params?.steps.containsKey("fragmentRecruitment") &&  params?.steps.fragmentRecruitment.containsKey("mashScreen"), \
	Channel.value([getModulePath(params.modules.fragmentRecruitment), "mash/paste"]),  pMashPasteChunk.out.sketch | collect(flat: false))

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
     PATH_IDX = 0
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

     ontReads | join(genomesMerged, by: SAMPLE_IDX) \
	| map { sample -> [sample[SAMPLE_IDX], sample[GENOMES_MERGED_IDX], sample[SAMPLE_2_IDX]] } \
	| set { minimapInput }

     pMinimap2(Channel.value(params?.steps.containsKey("fragmentRecruitment")), \
	Channel.value([Utils.getModulePath(params.modules.fragmentRecruitment), \
        "readMapping/minimap",  params.steps?.fragmentRecruitment?.mashScreen?.additionalParams?.minimap, \
        params.steps?.fragmentRecruitment?.mashScreen?.additionalParams?.samtoolsViewMinimap, \
	params.steps.containsKey("fragmentRecruitment")]), minimapInput)

      pBowtie2.out.mappedReads | mix(pBwa.out.mappedReads) | set { mappedReads }

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
     PATH_IDX = 0
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
	| join(Channel.value(DO_NOT_SET_IDENTITY_AUTOMATICALLY), by: SAMPLE_IDX) \
        | set { covermBowtieReadsInput  }

     minimapMappedReads \
	| join(covermGenomesInput, by: SAMPLE_IDX) \
	| join(Channel.value(ontMedianQuality), by: SAMPLE_IDX) \
        | set { covermMinimapReadsInput }

     covermBowtieReadsInput | mix(covermMinimapReadsInput) | pCovermCount

     // Found genomes reported by coverm are just file names. These lines join unique file names with file paths. 
     pCovermCount.out.foundGenomes |  splitCsv(sep: '\t', header: false) \
	| map { genomes -> genomes[GENOMES_IDX] } | flatten | join(genomesMap) | map { genome -> genome[GENOMES_IDX] } \
	| set { covermFilteredGenomes }

     covermFilteredGenomes | unique { genome -> genome.name } \
        | collect | map { genomes -> [ "EXTERNAL_GENOMES", genomes] } | set { foundGenomes }

     // Get Bin coverage statistics of the alignment
     covermFilteredGenomes | map { genome -> ["EXTERNAL_GENOMES", file(genome).name, genome] } | set { genomesSeperated }

     pCovermCount.out.foundGenomes | splitCsv(sep: '\t', header: false) \
        | map { genomes -> Utils.flattenTuple(genomes).flatten()} \
	| combine(genomesMap | map { genome -> genome.reverse() }, by: GENOMES_IDX) \
	| map { genomes -> [genomes[SAMPLE_2_IDX], genomes[GENOMES_2_IDX]]} \
	| groupTuple(by: SAMPLE_IDX) | set { foundGenomesInGroup }
  
     foundGenomesInGroup | pGenomeContigMapping

     pGenomeContigMapping.out.mapping | join(mappedShortReads, by: SAMPLE_IDX) \
        | combine(Channel.from("external")) | join(foundGenomesInGroup, by: SAMPLE_IDX) \
        | join(Channel.value(DO_NOT_SET_IDENTITY_AUTOMATICALLY), by: SAMPLE_IDX) \
        | set { shortReadMappingStatsInput }


     pGenomeContigMapping.out.mapping | join(minimapMappedReads, by: SAMPLE_IDX) \
        | combine(Channel.from("external")) | join(foundGenomesInGroup, by: SAMPLE_IDX) \
        | join(ontMedianQuality, by: SAMPLE_IDX)
        | set { minimapMappingStatsInput }

     pGetBinStatistics(getModulePath(params.modules.fragmentRecruitment), shortReadMappingStatsInput | mix(minimapMappingStatsInput))

  emit:
     genomesSeperated = genomesSeperated
     binsStats = pGetBinStatistics.out.binsStats     
     foundGenomes = foundGenomes
}

process pUnzip {

  label 'tiny'

  container "${params.ubuntu_image}"

  input:
  path x

  output:
  path("out/${x.baseName}${concatEnding}")

  script:
  ending = file(x).name.substring(file(x).name.lastIndexOf(".")) 
  concatEnding =  ending == ".gz" ? "" : ending 
  """
  mkdir out
  < $x zcat --force > out/${x.baseName}${concatEnding}
  """
}

workflow _wMashScreen {
   take: 
     sampleReads
     singleReads
     ontReads
     ontMedianQuality
   main:
     PATH_IDX = 0
     SAMPLE_IDX = 0
     BIN_ID_IDX = 1
     PATH_2_IDX = 2
     STATS_PATH = 1

     genomes = Channel.empty()
     if(params?.steps.containsKey("fragmentRecruitment") \
	&& params?.steps.fragmentRecruitment.containsKey("mashScreen") \
	&& params?.steps.fragmentRecruitment.mashScreen.containsKey("genomes")){
       Channel.fromPath(params?.steps.fragmentRecruitment?.mashScreen?.genomes) | splitCsv(sep: '\t', header: true) \
         | map { line -> file(line.PATH)} | pUnzip | set { genomes }
     }

     UNIQUE_IDX=0
     ALL_IDX=0
     genomes | map { genome -> genome.name } | unique | count | combine(genomes | count) \
	| filter { count -> count[UNIQUE_IDX]!=count[ALL_IDX] } | view { error "Genome file names must be unique!" } 

     genomes | map { genome -> [file(genome).name, genome] } | set { genomesMap }

     // Run mash and return matched genomes per sample
     _wRunMash(sampleReads | mix(ontReads), singleReads, genomesMap)

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
     genomes = _wGetStatistics.out.foundGenomes
     binsStats = binMap
     genomesSeperated = _wGetStatistics.out.genomesSeperated
}
