include { wSaveSettingsList } from '../config/module'

import java.lang.Math;


def getOutput(SAMPLE, RUNID, TOOL, filename){
    return SAMPLE + '/' + RUNID + '/' + params.modules.assembly.name + '/' +
          params.modules.assembly.version.major + "." +
          params.modules.assembly.version.minor + "." +
          params.modules.assembly.version.patch +
          '/' + TOOL + '/' + filename
}

def timestamp = new java.util.Date().format( 'YYYYMMdd-HHmmss-SSS')

/*
* This process uses kmer frequencies and the nonpareil diversity index to predict peak memory consumption on an assembler.
* A random forest regression model was trained on thousands of datasets of different environments in order to be able to predict
* the peak memory consumption.
*/
process pPredictFlavor {

    label 'tiny'

    tag "Sample: $sample"

    publishDir params.output, mode: "${params.publishDirMode}", saveAs: { filename -> getOutput("${sample}", params.runid, "predictFlavor", filename) }

    when params?.steps.containsKey("assembly") && params?.steps?.assembly.containsKey("megahit")

    container "${params.assemblerResourceEstimator_image}"

    input:
    val(modelType)
    tuple val(sample), path(interleavedReads), path(unpairedReads), val(nonpareilDiversity), path(kmerFrequencies71), path(kmerFrequencies21), path(kmerFrequencies13), path(model)
    val(minMemory)

    output:
    tuple file(".command.sh"), file(".command.out"), file(".command.err"), file(".command.log")
    tuple val("${sample}"), env(MEMORY), emit: memory
    tuple val("${sample}"), path("*.tsv"), emit: details

    shell:
    error = modelType == "sensitive" ?12 : 5
    '''
    zcat !{interleavedReads} !{unpairedReads} | seqkit stats --all -T > seqkit.stats.tsv
    GC_CONTENT=$(cut -d$'\t' -f 16 seqkit.stats.tsv | tail -n 1)
    MAX_READ_LEN=$(cut -d$'\t' -f 7 seqkit.stats.tsv | tail -n 1)
    if [ "${MAX_READ_LEN}" -gt 71 ]; then
    	MEMORY=$(cli.py predict -m !{model} -k21 !{kmerFrequencies21} -k71 !{kmerFrequencies71} -k13 !{kmerFrequencies13} -e !{error} -g ${GC_CONTENT} -d !{nonpareilDiversity} -o .)
    else
        MEMORY="!{minMemory}"
    fi
    echo -e "SAMPLE\tMEMORY" > !{sample}_ram_prediction.tsv
    echo -e "!{sample}\t${MEMORY}" >> !{sample}_ram_prediction.tsv
    '''
}


/*
*
* This method returns the maximum of the specified resources.
*/
def getMaxAvailableResource(type){
  return params.resources.subMap(params.resources.keySet()).values().collect({it[type]}).max()
}

/*
* 
* This closure returns the next higher cpu or memory value for specified exit codes of a failed tool run (e.g. due to memory restrictions).
* If the exit code is expected, this closure returns the next higher cpu/memory value based on the
* flavor with the next higher memory value.
* Initial CPU and memory values can be set by providing a predicted memory value. Based on the predicted memory value the next higher memory label
* of a user defined label is set as initial memory value. 
* The highest possible cpu/memory value is restricted by the highest cpu/memory value of all flavors. 
*/
def getNextHigherResource = { exitCodes, exitStatus, resourceType, attempt, memory, tool, defaults, sample ->  

                      // Check if Megahit failed based on a known exit code
                      if(exitStatus in exitCodes ){
     		         def currentResource = getResources(memory, tool, defaults, sample, attempt);
                         def currentResourceType = currentResource[resourceType];

                         return currentResourceType;
                      } else {
                         // if the exit code is not known then the memory value should not be increased
                         FLAVOR_INDEX_CONSTANT_TO_ADD = 1
                         def currentResource = getResources(memory, tool, defaults, sample, FLAVOR_INDEX_CONSTANT_TO_ADD);
                         def currentResourceType = currentResource[resourceType];

                         return currentResourceType;
                      }
}


process pMegahit {

    tag "Sample: $sample"

    cpus { getNextHigherResource([-9, 137, 247], task.exitStatus, "cpus", task.attempt, \
	"${memory}", params.steps.assembly.megahit, \
	params.resources.highmemLarge, "${sample}") }

    memory { getNextHigherResource([-9, 137, 247], task.exitStatus, "memory", task.attempt, \
	"${memory}", params.steps.assembly.megahit, \
	params.resources.highmemLarge, "${sample}") + ' GB' }

    publishDir params.output, mode: "${params.publishDirMode}", saveAs: { filename -> getOutput("${sample}", params.runid, "megahit", filename) }

    when params?.steps.containsKey("assembly") && params?.steps?.assembly.containsKey("megahit")

    container "${params.megahit_image}"

    input:
    tuple val(sample), path(interleavedReads, stageAs: 'interleaved.fq.gz'), path(unpairedReads), val(memory)

    output:
    tuple val("${sample}"), path("${sample}_contigs.fa.gz"), optional: true, emit: contigs
    tuple val("${sample}"), path("${sample}_contigs_stats.tsv"), optional: true, emit: contigsStats
    tuple val("${sample}"), path("${sample}_contigs.fastg"), env(maxKmer), optional: true, emit: fastg
    tuple file(".command.sh"), file(".command.out"), file(".command.err"), file(".command.log")

    shell:
    includeUnpairedReads = unpairedReads.name != "NOT_SET" ? " -r ${unpairedReads} " : ''
    convertToFastg = params?.steps?.assembly?.megahit?.fastg ? "TRUE" : "FALSE"
    template 'megahit.sh'
}


process pMetaspades {

    label 'highmemLarge'

    tag "$sample"

    publishDir params.output, mode: "${params.publishDirMode}", saveAs: { filename -> getOutput("${sample}", params.runid, "metaspades", filename) }

    memory { Utils.getMemoryResources(params.resources.highmemLarge, "${sample}", task.attempt, params.resources) }

    cpus { Utils.getCPUsResources(params.resources.highmemLarge, "${sample}", task.attempt, params.resources) }

    when params?.steps.containsKey("assembly") && params?.steps?.assembly.containsKey("metaspades")

    container "${params.metaspades_image}"

    input:
    tuple val(sample), path(interleavedReads, stageAs: 'interleaved.fq.gz')

    output:
    tuple val("${sample}"), path("${sample}_contigs.fa.gz"), emit: contigs
    tuple val("${sample}"), path("${sample}_contigs_stats.tsv"), emit: contigsStats
    tuple val("${sample}"), path("${sample}_contigs.fastg"), env(maxKmer), emit: fastg, optional: true
    tuple file(".command.sh"), file(".command.out"), file(".command.err"), file(".command.log")

    shell:
    outputFastg = params?.steps?.assembly?.metaspades?.fastg ? "TRUE" : "FALSE"
    memory = task.memory.toGiga()
    template 'metaspades.sh'
}


/*
*
* This workflow is only meant for testing purposes.
*
*/
workflow wTestMemSelection {
  main:
        exitCode = -9
        attempt = 2
        flavorRAMPredicted = "34"
        resourcesFlavor =  params.resources.highmemXLarge
        outCPUs = getNextHigherResource([-9, 137, 247], exitCode, "cpus", attempt, \
        flavorRAMPredicted, params.steps.assembly.megahit, \
        resourcesFlavor, "SAMPLE")

        outRAM = getNextHigherResource([-9, 137, 247], exitCode, "memory", attempt, \
        flavorRAMPredicted, params.steps.assembly.megahit, \
        resourcesFlavor, "SAMPLE")

        println("CPUs: " + outCPUs)
        println("RAM: " + outRAM)
}



/*
 * Input: 
 *  - List of the format [SAMPLE, Reads Paired, Reads Unpaired] 
 *  - Channel containing the values [SAMPLE, nonpareil tsv file]
 *  - Channel containing kmer frequency file and sample identifier [SAMPLE, kmer frequency]
 *
 * Output is of the format [SAMPLE, CONTIGS] for contigs and [SAMPLE, fastg] for fastq files.
 * 
 */
workflow wShortReadAssemblyList {
    take:
       readsList
       nonpareil
       kmerFrequencies
    main:
       _wAssembly(readsList, nonpareil, kmerFrequencies)
    emit:
      contigs = _wAssembly.out.contigs
      fastg = _wAssembly.out.fastg
}



/*
 * Takes two tab separated file of files containing paired and optional single reads 
 * as input and produces assembly results.
 * Input files must have two columns seperated by tabs:
 * SAMPLE and READS
 *
 * Output is of the format [SAMPLE, CONTIGS]
 * 
 */
workflow wShortReadAssemblyFile {
    main:
       SAMPLE_IDX = 0       
       SAMPLE_PAIRED_IDX = 1
       UNPAIRED_IDX = 2

       readsPaired = Channel.empty()
       if(params.steps.assembly.input.containsKey("paired")) {
       	 Channel.from(file(params.steps.assembly.input.paired)) | splitCsv(sep: '\t', header: true) \
             | map { it -> [ it.SAMPLE, it.READS]} | set { readsPaired  }
       }

       readsSingle = Channel.empty()
       if(params.steps.assembly.input.containsKey("single")) {
         Channel.from(file(params.steps.assembly.input.single)) | splitCsv(sep: '\t', header: true) \
             | map { it -> [ it.SAMPLE, it.READS]} | set { readsSingle  }
       }

       readsPaired | join(readsSingle, by: SAMPLE_IDX, remainder: true) \
	| map { sample -> sample[UNPAIRED_IDX] == null ? \
		[sample[SAMPLE_IDX], sample[SAMPLE_PAIRED_IDX], file("NOT_SET")] : sample } \
	| set { reads }

        wSaveSettingsList(reads | map { it -> it[SAMPLE_IDX] })

       _wAssembly(reads, Channel.empty(), Channel.empty())
    emit:
      contigs = _wAssembly.out.contigs
}

/*
* This method sets minimum memory and cpu if the predicted memory is below a user defined threshold.
*
*/
def setMinLabel(assembler, labelMap, memoryLabelMap, sortedMemorySet, nextHigherMemoryIndex){
       def label = memoryLabelMap[sortedMemorySet[nextHigherMemoryIndex]]

       def cpus =  labelMap[label]["cpus"]
       def memory =  labelMap[label]["memory"]
       def minLabel = assembler.resources.RAM.predictMinLabel

       if(minLabel in labelMap.keySet()) {
           def minLabelCpus = labelMap[minLabel]["cpus"]
           def minLabelMemory = labelMap[minLabel]["memory"]
            if(minLabelCpus > cpus || minLabelMemory > memory){
	        return ["cpus": minLabelCpus, "memory": minLabelMemory];
	    } else {
	        return ["cpus": cpus, "memory": memory];
            }
       } else {
            println("WARNING: Unknown predictMinLabel parameter specified! Setting next higher resource flavor");
	    return ["cpus": cpus, "memory": memory];
       }
}


/*
*
* This method sets memory and ram for an assembler based on 
* user defined parameter ('PREDICT' or 'DEFAULT'). 
*
*/
def getResources(predictedMemory, assembler, defaults, sample, attempt){

     // get map of memory values and flavor names (e.g. key: 256, value: large)
     def memoryLabelMap = params.resources.findAll().collectEntries( { [it.value.memory, it.key] });

     // Get map of resources with label as key
     def labelMap = params.resources.findAll()\
		.collectEntries( { [it.key, ["memory" : it.value.memory, "cpus" : it.value.cpus ] ] })

     // get memory values as list
     def memorySet = memoryLabelMap.keySet().sort()



     switch(assembler.resources.RAM.mode){
       case 'PREDICT':
     	  def predictedMemoryCeil = (int)Math.ceil(Float.parseFloat(predictedMemory))

          // add predicted memory to list and sort to get next index of next higher resource label memory
     	  def updatedMemorySet = memorySet + predictedMemoryCeil
     	  def sortedUpdatedMemorySet = updatedMemorySet.sort()
     	  def predictedMemoryIndex = sortedUpdatedMemorySet.findIndexOf({ it == predictedMemoryCeil })
          
     	  def nextHigherMemoryIndex = predictedMemoryIndex + attempt

     	  if(nextHigherMemoryIndex  >= sortedUpdatedMemorySet.size()){
             // In case it is already the highest possible memory setting
             // then try the label with the highest memory
             println("Warning: Predicted Memory " + predictedMemoryCeil \
		+ " of the dataset " + sample + " is greater or equal to the flavour with the largest RAM specification.")

             def label = memoryLabelMap[memorySet[memorySet.size() -1]];
             def cpus =  labelMap[label]["cpus"];
             def ram =  labelMap[label]["memory"];

	     return ["cpus": cpus, "memory": ram];
	  } else {
               // If the memory value is not the highest possible then check if the the user provided 
               // flavor with the minimum memory value should be provided
	       switch(assembler.resources.RAM.predictMinLabel) {
        	  case "AUTO":
       			def label = memoryLabelMap[sortedUpdatedMemorySet[nextHigherMemoryIndex]];
		        def cpus = labelMap[label]["cpus"];
		        def memory = labelMap[label]["memory"];
			return ["cpus": cpus, "memory": memory];
         	  default:
	    		def minLabel = setMinLabel(assembler, labelMap, memoryLabelMap, sortedUpdatedMemorySet, nextHigherMemoryIndex);
                        return minLabel
       	       }
       	  }
          break;
       case 'DEFAULT':

          // If the memory value is not predicted do the following
          def defaultCpus = defaults.cpus

          def defaultMemory = defaults.memory

     	  def defaultMemoryIndex = memorySet.findIndexOf({ it == defaultMemory })

 	  def nextHigherMemoryIndex = defaultMemoryIndex + attempt - 1 


     	  if(nextHigherMemoryIndex >= memorySet.size()){
             // In case it is already the highest possible memory setting
             // then try the label with the highest memory
             println("Warning: Highest possible memory setting " \
		+ " of the dataset " + sample + " is already reached." \
                + " Retry will be done using a flavor with the highest possible RAM configuration")

             def label = memoryLabelMap[memorySet[memorySet.size()-1]]
             def cpus =  labelMap[label]["cpus"]
             def ram =  labelMap[label]["memory"]

	     return ["cpus": cpus, "memory": ram];
	  } else {
       	     def label = memoryLabelMap[memorySet[nextHigherMemoryIndex]];
             def cpus = labelMap[label]["cpus"];
             def memory = labelMap[label]["memory"];
	     return ["cpus": cpus, "memory": memory];
       	  }

          return ["cpus": cpus, "memory": memory] 
     }
}


/*
* Get memory and cpus of the minimum specified resource flavor
*/
def getMinResources(){ 
        minLabel = params.steps.assembly.megahit.resources.RAM.predictMinLabel
        minCPUs = params.resources[minLabel].cpus
        minMem = params.resources[minLabel].memory
        return ["memory": minMem, "cpus": minCPUs]
}

workflow _wCalculateMegahitResources {
       take:
         readsList
         nonpareil
         kmerFrequencies
       main:
         SAMPLE_IDX = 0
         NONPAREIL_METRICS_IDX = 1
 
         // figure out whether resources should be predicted or not
         readsList | branch {
        	predict: params?.steps?.assembly?.megahit?.resources?.RAM?.mode == "PREDICT"
        	doNotPredict: params?.steps?.assembly?.megahit?.resources?.RAM?.mode == "DEFAULT"
         } | set { resourceType }

         model = Channel.empty()
         modelType = Channel.empty()
         if(params.steps.containsKey("assembly") && params.steps.assembly.containsKey("megahit") \
		&& params?.steps?.assembly?.megahit?.additionalParams.contains("meta-sensitive")){
         	model = Channel.value(file("${projectDir}/models/assembler/megahit/sensitive.pkl"))
		modelType = Channel.value("sensitive")
	 } else {
		modelType = Channel.value("default")
         	model = Channel.value(file("${projectDir}/models/assembler/megahit/default.pkl"))
	 }
   
         resourceType.predict | join(nonpareil | splitCsv(header: true, sep: '\t') \
	  | map{ it -> [it[SAMPLE_IDX], it[NONPAREIL_METRICS_IDX].diversity]  }) \
          | join(kmerFrequencies) | combine(model) | map { dataset -> dataset.flatten() } | set { predictFlavorInput }
         
         // Get memory of the minimum specified resource label
         minResources = getMinResources()
         pPredictFlavor(modelType, predictFlavorInput, Channel.value(minResources.memory))

         PREDICTED_RAM_IDX = 1

         pPredictFlavor.out.memory \
          | collectFile(newLine: true, seed: "SAMPLE\tPREDICTED_RAM", storeDir: params.logDir){ item ->
        	[ "predictedMegahitRAM." + timestamp + ".tsv", item[SAMPLE_IDX] + '\t' + item[PREDICTED_RAM_IDX]  ]
    	  }

         resourceType.doNotPredict | map{ it -> it + "NoPrediction" } \
	   | mix(readsList | join(pPredictFlavor.out.memory)) | set { resources }
       emit:
         resources
}


/*
*
* Input:
*  - List of the format [SAMPLE, Reads Paired, Reads Unpaired] 
*  - Channel containing the values [SAMPLE, nonpareil tsv file]
*  - Channel containing kmer frequency file and sample identifier [SAMPLE, kmer frequency]
*/
workflow _wAssembly {
     take:
       readsList
       nonpareil
       kmerFrequencies
     main:
       _wCalculateMegahitResources(readsList, nonpareil, kmerFrequencies) | pMegahit

       // Metaspades does only accept paired end
       // Thats why orphaned reads are filtered out
       PAIRED_END_ID = 1
       SAMPLE_IDX = 0
       readsList | map { seq -> [seq[SAMPLE_IDX], seq[PAIRED_END_ID]] } | pMetaspades

       pMegahit.out.contigs | mix(pMetaspades.out.contigs) | set { contigs }

       pMegahit.out.fastg | mix(pMetaspades.out.fastg) | set { fastg }
       
    emit:
      contigs = contigs
      fastg = fastg
}
