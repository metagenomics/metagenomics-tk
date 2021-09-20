nextflow.enable.dsl=2

def getOutput(RUNID, TOOL, filename){
    return "AGGREGATED" + '/' +  RUNID + '/' + params.modules.cooccurrence.name + '/' + 
         params.modules.cooccurrence.version.major + "." +  
         params.modules.cooccurrence.version.minor + "." +  
         params.modules.cooccurrence.version.patch +  
         '/' + TOOL + '/' + filename
}


process pVerticalConcat {

    errorStrategy 'ignore'

    label 'tiny'

    when params.steps.containsKey("cooccurrence")

    input:
    file('sample?')

    output:
    path("abundance.tsv")

    shell:
    template("verticalConcat.sh")
}


process pBuildNetwork {

    errorStrategy 'ignore'

    label 'large'
   
    container "pbelmann/cooccurrence:${params.cooccurrence_tag}"

    publishDir params.output, saveAs: { filename -> getOutput(params.runid, "matrix", filename) }

    when params.steps.containsKey("cooccurrence")

    input:
    file('abundance.tsv')
    file('gtdb.tsv')

    output:
    path("community.tsv"), emit: community
    path("output.graphml"), emit: graphml

    shell:
    template("coocurrence.sh")
}


process pVerticalConcatFinal {

    errorStrategy 'ignore'

    label 'tiny'

    publishDir params.output, saveAs: { filename -> getOutput(params.runid, "matrix", filename) }

    when params.steps.containsKey("cooccurrence")

    input:
    file('sample?')

    output:
    path("abundance.tsv")

    shell:
    template("verticalConcat.sh")
}


/*
*
* This workflow builds cooccurrence network based on abundance information and gtdb taxonomy assignment.
* Input:
*   * Count: Abundance matrix file containing MAG abundance values per sample
*   * GTDB: GTDB assignments per sample
*
*/
workflow wCooccurrenceFile {
  main: 
    Channel.from(file(params.steps.cooccurrence.count))\
       | splitCsv(sep: '\t', header: false, skip: 1) \
       | set {count}
    _wCooccurrence(count, Channel.from(file(params.steps.cooccurrence.gtdb)))
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
*
*/
workflow wCooccurrenceList {
  take:
    count
    gtdb
  main:
     gtdb | collectFile(keepHeader: true) { file, dataset -> ["gtdb", file.text]}\
       | set {gtdbConcatenated}
    _wCooccurrence(count, gtdbConcatenated)
}


workflow _wCooccurrence {
   take: 
     count
     gtdbConcatenated
   main:
     COLLECT_BUFFER=10000

     // Collect abundance/ count information of every sample
     count | map { id, count -> file(count) }\
        | buffer( size: COLLECT_BUFFER, remainder: true ) \
        | pVerticalConcat | collect \
        | pVerticalConcatFinal | set { abundance }

     // Run the network process
     pBuildNetwork(abundance, gtdbConcatenated)
}
