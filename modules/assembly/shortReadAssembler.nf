nextflow.enable.dsl=2

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
    tuple val(sample), path(interleavedReads), path(unpairedReads), val(nonpareilDiversity), path(kmerFrequencies), path(model)

    output:
    tuple file(".command.sh"), file(".command.out"), file(".command.err"), file(".command.log")
    tuple val("${sample}"), env(MEMORY), emit: memory

    shell:
    '''
    READS_COUNTER=$(zcat !{interleavedReads} !{unpairedReads} | awk '{s++}END{print s/4}')
    BASEPAIRS_COUNTER=$(zcat !{interleavedReads} !{unpairedReads} | seqkit stats -T | cut -d$'\t' -f 5 | tail -n 1)
    MEMORY=$(cli.py predict -m !{model} -k !{kmerFrequencies} -b ${BASEPAIRS_COUNTER} -d !{nonpareilDiversity} -r ${READS_COUNTER})
    echo -e "Memory: ${MEMORY}\nBasepairs: ${BASEPAIRS_COUNTER}\nReads: ${READS_COUNTER}" 
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
* This closue returns the next higher cpu or memory value for specified exit codes of a failed tool run (e.g. due to memory restrictions).
* If the exit code is expected, this closure returns the next higher cpu/memory value computed by 
* the formula 2^(number of attempts) * (cpu/memory value of the assigned flavour).
* Initial CPU and memory values can be set by providing a predicted memory value. Based on the predicted memory value the next higher memory label
* of a user defined label is set as initial memory value. 
* The highest possible cpu/memory value is restricted by the highest cpu/memory value of all flavors. 
*/
def getNextHigherResource = { exitCodes, exitStatus, resourceType, attempt, memory, tool, defaults, sample ->  

                      maxResource = getMaxAvailableResource(resourceType);

		      currentResource = getResources(memory, tool, defaults, sample)

                      currentResourceType = currentResource[resourceType]

                      if(exitStatus in exitCodes ){

                        if(currentResourceType * attempt < maxResource){
                          return Math.pow(2, attempt - 1) * currentResourceType as long;
                        } else {
                          return maxResource;
                        }
                      } else {
                          return currentResourceType;
                      }
}


process pMegahit {

    tag "Sample: $sample"

    cpus { getNextHigherResource([-9, 137, 247], task.exitStatus, "cpus", task.attempt, \
	"${memory}", params.steps.assembly.megahit, \
	params.resources.large, "${sample}") }

    memory { getNextHigherResource([-9, 137, 247], task.exitStatus, "memory", task.attempt, \
	"${memory}", params.steps.assembly.megahit, \
	params.resources.large, "${sample}") + ' GB' }

    publishDir params.output, mode: "${params.publishDirMode}", saveAs: { filename -> getOutput("${sample}", params.runid, "megahit", filename) }

    when params?.steps.containsKey("assembly") && params?.steps?.assembly.containsKey("megahit")

    container "${params.megahit_image}"

    input:
    tuple val(sample), path(interleavedReads, stageAs: 'interleaved.fq.gz'), path(unpairedReads), val(memory)

    output:
    tuple val("${sample}"), path("${sample}_contigs.fa.gz"), emit: contigs
    tuple val("${sample}"), path("${sample}_contigs_stats.tsv"), emit: contigsStats
    tuple val("${sample}"), path("${sample}_contigs.fastg"), env(maxKmer), emit: fastg, optional: true
    tuple file(".command.sh"), file(".command.out"), file(".command.err"), file(".command.log")

    shell:
    includeUnpairedReads = unpairedReads.name != "NOT_SET" ? " -r ${unpairedReads} " : ''
    convertToFastg = params?.steps?.assembly?.megahit?.fastg ? "TRUE" : "FALSE"
    template 'megahit.sh'
}


process pMetaspades {

    label 'large'

    tag "$sample"

    publishDir params.output, mode: "${params.publishDirMode}", saveAs: { filename -> getOutput("${sample}", params.runid, "metaspades", filename) }

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
    template 'metaspades.sh'
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

       _wAssembly(reads, Channel.empty(), Channel.empty())
    emit:
      contigs = _wAssembly.out.contigs
}

/*
* This method sets minimum memory and cpu if the predicted memory is below a user defined threshold.
*
*/
def setMinLabel(assembler, labelMap, memoryLabelMap, sortedMemorySet, nextHigherMemoryIndex){
       label = memoryLabelMap[sortedMemorySet[nextHigherMemoryIndex]]
       cpus =  labelMap[label]["cpus"]
       memory =  labelMap[label]["memory"]
       minLabel = assembler.resources.RAM.predictMinLabel

       if(minLabel in labelMap.keySet()) {
           minLabelCpus = labelMap[minLabel]["cpus"]
           minLabelMemory = labelMap[minLabel]["memory"]
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
def getResources(predictedMemory, assembler, defaults, sample){
     switch(assembler.resources.RAM.mode){
       case 'PREDICT':
          // get map of memory values and flavor names (e.g. key: 256, value: large)
     	  memoryLabelMap = params.resources.findAll().collectEntries( { [it.value.memory, it.key] })

          // Get map of resources with label as key  
     	  labelMap = params.resources.findAll()\
		.collectEntries( { [it.key, ["memory" : it.value.memory, "cpus" : it.value.cpus ] ] })

          // get memory values as list
     	  memorySet = memoryLabelMap.keySet()

     	  predictedMemoryCeil = Math.ceil(Float.parseFloat(predictedMemory))

          // add predicted memory to list and sort to get next index of next higher resource label memory 
     	  memorySet = memorySet + predictedMemoryCeil
     	  sortedMemorySet = memorySet.sort()
     	  predictedMemoryIndex = sortedMemorySet.findIndexOf({ it == predictedMemoryCeil })

     	  nextHigherMemoryIndex = predictedMemoryIndex + 1
     	  if(nextHigherMemoryIndex  == sortedMemorySet.size()){
             // In case it is already the highest possible memory setting
             // then try the label with the highest memory
             println("Warning: Predicted Memory " + predictedMemoryCeil \
		+ " of the dataset " + sample + " is greater or equal to the flavour with the largest RAM specification.")
             label = memoryLabelMap[sortedMemory[predictedMemoryIndex]]
             cpus =  labelMap[label]["cpus"]
	     return ["cpus": cpus, "memory": sortedMemorySet[predictedMemoryIndex]];
	  } else {
	       switch(assembler.resources.RAM.predictMinLabel) {
        	  case "AUTO":
       			label = memoryLabelMap[sortedMemorySet[nextHigherMemoryIndex]]
		        cpus = labelMap[label]["cpus"]
		        memory = labelMap[label]["memory"]
			return ["cpus": cpus, "memory": memory];
         	  default:
	    		return setMinLabel(assembler, labelMap, memoryLabelMap, sortedMemorySet, nextHigherMemoryIndex);
       	       }
       	  }
          break;
       case 'DEFAULT':
          cpus = defaults.cpus
          memory = defaults.memory
          return ["cpus": cpus, "memory": memory] 
     }
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
         if(params.steps.containsKey("assembly") && params.steps.assembly.containsKey("megahit") \
		&& params?.steps?.assembly?.megahit?.additionalParams.contains("meta-sensitive")){
         	model = Channel.value(file("${baseDir}/models/assembler/megahit/sensitive.pkl"))
	 } else {
         	model = Channel.value(file("${baseDir}/models/assembler/megahit/default.pkl"))
	 }
   
         resourceType.predict | join(nonpareil | splitCsv(header: true, sep: '\t') \
	  | map{ it -> [it[SAMPLE_IDX], it[NONPAREIL_METRICS_IDX].diversity]  }) \
          | join(kmerFrequencies) | combine(model) | pPredictFlavor

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

       if(params.summary){
         pMegahit.out.contigsStats | mix(pMetaspades.out.contigsStats) | collectFile(newLine: false, keepHeader: true, storeDir: params.output + "/summary/" ){ item ->
           [ "contigs_stats.tsv", item[1].text ]
         }
       }
      
       pMegahit.out.contigs | mix(pMetaspades.out.contigs) | set { contigs }

       pMegahit.out.fastg | mix(pMetaspades.out.fastg) | set { fastg }
       
    emit:
      contigs = contigs
      fastg = fastg
}
