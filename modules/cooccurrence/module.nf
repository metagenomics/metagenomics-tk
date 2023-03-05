nextflow.enable.dsl=2

include { pDumpLogs } from '../utils/processes'

def getOutput(RUNID, TOOL, filename){
    return "AGGREGATED" + '/' +  RUNID + '/' + params.modules.cooccurrence.name + '/' + 
         params.modules.cooccurrence.version.major + "." +  
         params.modules.cooccurrence.version.minor + "." +  
         params.modules.cooccurrence.version.patch +  
         '/' + TOOL + '/' + filename
}


process pVerticalConcat {

    label 'tiny'

    cache 'deep'

    when params.steps.containsKey("cooccurrence")

    input:
    file('sample?')

    output:
    path("abundance.tsv")

    shell:
    template("verticalConcat.sh")
}


process pBuildNetwork {

    label 'large'

    cache 'deep'
   
    container "${params.cooccurrence_image}"

    publishDir params.output, mode: "${params.publishDirMode}", saveAs: { filename -> getOutput(params.runid, "network", filename) }

    when params.steps.containsKey("cooccurrence") && params.steps.cooccurrence.containsKey("inference")

    input:
    file('abundance.tsv')
    file('gtdb.input.tsv')

    output:
    path("community.tsv"), emit: community
    path("edges_index.tsv"), emit: edges
    path("gtdb.tsv"), emit: gtdb
    path("output.graphml"), emit: graph
    path("output_raw.graphml"), emit: graphRaw
    tuple file(".command.sh"), file(".command.out"), file(".command.err"), file(".command.log")

    shell:
    batchSize = params.steps.cooccurrence.containsKey("metabolicAnnotation") ? params.steps.cooccurrence.metabolicAnnotation.additionalParams.metabolicEdgeBatches : 1
    template("coocurrenceCreate.sh")
}


process pUpdateNetwork {

    label 'tiny'

    container "${params.cooccurrence_image}"

    publishDir params.output, mode: "${params.publishDirMode}",\
	saveAs: { filename -> getOutput(params.runid, "network", filename) }

    when params.steps.containsKey("cooccurrence") && params.steps.cooccurrence.containsKey("metabolicAnnotation")

    input:
    path(graph)
    path(attributes)

    output:
    path("updated.graphml"), emit: updatedGraphml
    tuple file(".command.sh"), file(".command.out"), file(".command.err"), file(".command.log")

    shell:
    """
    Rscript /cooccurrence.R update -g !{graph} -e !{attributes} -o .
    """
}

process pVerticalConcatFinal {

    label 'tiny'

    cache 'deep'

    publishDir params.output, mode: "${params.publishDirMode}", saveAs: { filename -> getOutput(params.runid, "matrix", filename) }

    when params.steps.containsKey("cooccurrence")

    input:
    file('sample?')

    output:
    path("abundance.tsv")

    shell:
    template("verticalConcat.sh")
}


MAX_SMETANA_RETRIES = 10
process pSmetanaEdges {

    label 'large'

    cache 'deep'

    tag "Edges Batch: ${edges}, Repeat: ${repeat}"

    errorStrategy { if(task.attempt <= MAX_SMETANA_RETRIES){ sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' } else { return 'ignore'  } }

    container "${params.smetana_image}"

    beforeScript params?.steps.containsKey("metabolomics") \
	? Utils.getBeforeScript(params?.steps?.cooccurrence?.beforeProcessScript.trim(), params.smetana_image) \
	: ""

    when params.steps.containsKey("cooccurrence") && params.steps.cooccurrence.containsKey("metabolicAnnotation")

    input:
    each repeat
    tuple path(edges), path(xmls)  

    output:
    path("edge_attributes.tsv"), emit: edgeAttributes
    tuple val("${edges}_${repeat}"), val(getOutput(params.runid, "smetana", "")), val(params.LOG_LEVELS.ALL), file(".command.sh"), \
	file(".command.out"), file(".command.err"), file(".command.log"), emit: logs

    shell:
    template "smetanaEdges.sh"
}


/*
*
* This workflow builds cooccurrence network based on abundance information and gtdb taxonomy assignment.
* Input:
*   * Count: Abundance matrix file containing MAG abundance values per sample
*   * GTDB: GTDB assignments per sample
*   * Models: Genome scale metabolic models of all MAGs 
*/
workflow wCooccurrenceFile {
  main: 
    Channel.from(file(params.steps.cooccurrence.input.count))\
       | splitCsv(sep: '\t', header: false, skip: 1) \
       | set {count}

    models = Channel.empty()
    if(params.steps.cooccurrence.input.containsKey("models")){
      Channel.from(file(params.steps.cooccurrence.input.models))\
        | splitCsv(sep: '\t', header: true) \
        | map { bin -> [bin.BIN_ID, file(bin.PATH)] } | set {models}
    }
    _wCooccurrence(count, Channel.from(file(params.steps.cooccurrence.input.gtdb)), models)

}


/*
* This helper method removes the suffix of file. 
* Example: test1_bin1.fa -> test1_bin1
*/
def fixModelID(modelName){
    if(modelName.contains(".")){
       return modelName.take(modelName.lastIndexOf('.'))
    } else {
       return modelName
    }
}

/*
*
* This workflow builds cooccurrence network based on abundance information and gtdb taxonomy assignment.
*
* Input:
*  * The count channel must contain entries of the format: [sampleid,  /path/to/abundance/file]
*    Abundance file must contain the following to columns:
*    Genome  SAMPLE_NAME
*    test1_bin.1     61.92912
*    test2_bin.2     115.32025
*  
*  * The gtdb channel contains gtdb assignments of every sample.
*   
*  * The models channel takes the model name and its file path as input:
*    [ERR2592253_bin.23, /meta_test/medium/cooccurrence/ERR2592253_bin.23.fa.xml]
*
*/
workflow wCooccurrenceList {
  take:
    count
    gtdb
    models
  main:
     gtdb | collectFile(keepHeader: true) { file, dataset -> ["gtdb", file.text]}\
       | set {gtdbConcatenated}

    MODEL_NAME_IDX = 0
    MODEL_PATH_IDX = 1
    models | map { model -> [fixModelID(model[MODEL_NAME_IDX]), model[MODEL_PATH_IDX]]} | set { fixedModel }
    _wCooccurrence(count, gtdbConcatenated, fixedModel)
}


workflow _wCooccurrence {
   take: 
     countCh
     gtdbConcatenated
     models
   main:
     COLLECT_BUFFER=10000

     // Collect abundance/ count information of every sample
     countCh | combine(countCh | count )  \
	| filter({ sample, abundance, count -> count > 1 }) \
	| map { sample, abundance, count -> file(abundance) }\
        | buffer( size: COLLECT_BUFFER, remainder: true ) \
        | pVerticalConcat | collect \
        | pVerticalConcatFinal | set { abundance }

     // Run the network process
     pBuildNetwork(abundance, gtdbConcatenated)

     MODEL_IDX = 0
     MODEL_1_PATH_IDX = 1
     MODEL_PATH_IDX = 3
     MODEL_2_IDX = 1
     BATCH_IDX = 2
     models | map { model -> model[MODEL_1_PATH_IDX] } \
	| collect(flat:false) \
	| collect(flat:false) \
	| set { modelList }

     // Get all edges, map node ids to model paths (e.g. ERR2592252_bin.25 to 
     // /meta_test/medium/cooccurrence/ERR2592252_bin.25.fa.model.xml)
     pBuildNetwork.out.edges |  splitCsv(header: true, sep: '\t') | map{ bin -> [bin.V1, bin.V2, bin.IDX]} \
	| combine(models, by: MODEL_IDX) \
        | map { bin -> [bin[MODEL_2_IDX], bin[MODEL_PATH_IDX], bin[BATCH_IDX]]} \
	| combine(models, by: MODEL_IDX) \
        | map { bins -> [bins[MODEL_1_PATH_IDX].name, bins[MODEL_PATH_IDX].name, bins[BATCH_IDX]]} \
        | set { rowEdges }

     // Create batches of edges in separate files
     MODEL_FILES = 2
     BATCH_NAME_IDX = 1
     rowEdges | collectFile(){ edge -> [ "edges_" + edge[BATCH_IDX], edge.take(MODEL_FILES).join('\t') + '\n' ] } \
	| map { f -> [f.name.split('_')[BATCH_NAME_IDX], f] } | set{batchedFileEdges}

     // Create Maps between Model Name, Model Path and BinId
     MODEL_NAME_IDX = 0
     models | map { model -> [model[MODEL_1_PATH_IDX].name, model[MODEL_NAME_IDX]] } | set { modelNameToBinIdMap }
     models | map { model -> [model[MODEL_1_PATH_IDX].name, model[MODEL_1_PATH_IDX]] } | set { modelNameToModelPathMap }

     // Models of batched edges are created
     // Example: [15, [/meta_test/medium/cooccurrence/ERR2592253_bin.11.fa.xml,  
     // /meta_test/medium/cooccurrence/ERR2592268_bin.45.fa.xml]]
     // 15 is the batch id
     MODEL_NAME_1_IDX = 1
     BATCH_NAME_1_IDX = 0
     MODEL_PATHS = 3
     rowEdges | combine(modelNameToModelPathMap, by: MODEL_NAME_IDX) | combine(modelNameToModelPathMap \
	| map { model -> model.reverse()}, by: MODEL_NAME_1_IDX) \
	| map { model -> model.takeRight(MODEL_PATHS) } | groupTuple(by: BATCH_NAME_1_IDX) \
	| map { edges -> [edges[BATCH_NAME_1_IDX], edges.tail().flatten().unique()] } \
	| set { batchedEdges }

     //Run Smetana
     pSmetanaEdges(Channel.from(1..params.steps.cooccurrence.metabolicAnnotation.additionalParams.metabolicEdgeReplicates),\
	batchedFileEdges | join(batchedEdges) | map{ model -> model.tail()})

     // Get edge metrics and update the previously creatd network
     METRICS = 6
     pSmetanaEdges.out.edgeAttributes | splitCsv(sep: '\t', skip: 1) \
	| combine(modelNameToBinIdMap, by: MODEL_NAME_IDX) \
	| combine(modelNameToBinIdMap | map { model -> model.reverse()}, by: MODEL_NAME_1_IDX) \
	| map { edge -> edge.takeRight(METRICS) } \
	| collectFile(storeDir: params.output + "/" + getOutput(params.runid, "smetana", ""), \
	seed: ["medium", "size", "mip", "mro", "V1", "V2"].join('\t'), \
	newLine: true){ line -> ['edgeAttributes.tsv', line.join("\t")] } \
	| set { edgeAttributes } 

     pSmetanaEdges.out.logs | pDumpLogs

     pUpdateNetwork(pBuildNetwork.out.graph, edgeAttributes)
}
