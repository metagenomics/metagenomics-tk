nextflow.enable.dsl=2


def getOutput(RUNID, TOOL, filename){
    return "AGGREGATED" + '/' +  RUNID + '/' + params.modules.dereplication.name + '/' + 
         params.modules.dereplication.version.major + "." +  
         params.modules.dereplication.version.minor + "." +  
         params.modules.dereplication.version.patch +  
         '/' + TOOL + '/' + filename
}


process pMashSketchGenome {

    container "${params.mash_image}"

    errorStrategy 'retry'

    label 'tiny'

    when params?.steps.containsKey("dereplication") &&  params?.steps.dereplication.containsKey("pasolli")

    input:
    tuple path("g.fa"), val(binid)

    output:
    path("${binid}.msh"), emit: sketch
    tuple env(GENOME_PATH), val("${binid}"), emit: stagedGenome

    shell:
    '''
    ln -s g.fa !{binid}
    mash sketch !{params.steps.dereplication.pasolli.additionalParams.mash_sketch} !{binid} -o !{binid}.msh
    GENOME_PATH=$(readlink -f g.fa)
    '''
}


process pMashPaste {

    container "${params.mash_image}"

    errorStrategy 'retry'

    label 'large'

    input:
    path sketches, stageAs: 'sketch*.msh'

    output:
    file('final_sketch.msh')

    shell:
    '''
    mash paste final_sketch !{sketches}
    '''
}


process pMashDist {

    container "${params.mash_image}"

    errorStrategy 'retry'

    label 'large'

    when params?.steps.containsKey("dereplication") &&  params?.steps.dereplication.containsKey("pasolli")

    input:
    path sketches, stageAs: 'sketch*.msh'

    output:
    file('distances.tsv')

    shell:
    '''
    mash paste final_sketch !{sketches}
    mash dist !{params.steps.dereplication.pasolli.additionalParams.mash_dist} -p !{task.cpus} final_sketch.msh final_sketch.msh | cut -f 1,2,3 > distances.tsv
    '''
}


process pClusterDistances {

    errorStrategy 'retry'

    input:
    file('distances.tsv')

    container "${params.python_env_image}"

    label 'medium'

    output:
    tuple file("distances.tsv"), file('out/clusters.tsv') 

    shell:
    '''
    mkdir out
    cluster.py -i distances.tsv !{params.steps.dereplication.pasolli.additionalParams.cluster} -o out
    '''
}


process pSelectRepresentative {

    errorStrategy 'retry'

    input:
    path genome_table
    tuple file("distance"), file("cluster")

    container "${params.python_env_image}"

    publishDir params.output, saveAs: { filename -> getOutput(params.runid, "pasolli/selectedRepresentatives", filename) }

    label 'medium'

    output:
    path("intermediate_clusters.tsv"), emit: clusters
    path("refinement/representatives_to_compare.tsv"), emit: representatives
    tuple file(".command.sh"), file(".command.out"), file(".command.err"), file(".command.log")

    shell:
    template 'selectRepresentative.sh'
}


process pANIb {

    errorStrategy 'retry'

    label 'small'

    input:
    file("genome1*") 
    file("genome2*") 

    container "${params.ani_image}"

    when:
    params.steps.dereplication.pasolli.method.contains("ANI")

    output:
    file("*.out/out.tsv") 

    shell:
    template 'anib.sh'
}


process pTETRA {

    label 'small'

    input:
    file("genome1*") 
    file("genome2*") 

    container "${params.ani_image}"

    output:
    file("*.out/out.tsv") 

    when:
    params.steps.dereplication.pasolli.method.contains("TETRA")

    shell:
    template 'tetra.sh'
}

process pGetCluster {

    label 'tiny'

    publishDir params.output, saveAs: { filename -> getOutput(params.runid, "pasolli/clusters", filename) }

    container "${params.python_env_image}"

    input:
    path cluster, stageAs: 'cluster'
    path ani_values, stageAs: 'ani_values'
    path genomeAttributes, stageAs: 'genomeAttributes'

    output:
    path 'final_clusters.tsv', emit: final_clusters
    path 'ani_values.tsv', emit: ani_values
    tuple file(".command.sh"), file(".command.out"), file(".command.err"), file(".command.log")

    shell:
    '''
    mkdir out
    cat <(echo "GENOME_A\tGENOME_B\tANI")  <(sed 's/ /\t/g' !{ani_values}) > ani_values.tsv
    get_final_cluster.py -i !{genomeAttributes} -c !{cluster} -r ani_values.tsv -o out  -a !{params.steps.dereplication.pasolli.additionalParams.representativeAniCutoff}
    cp out/representatives.tsv final_clusters.tsv
    '''
}

process pFinalize {

    input:
    val finalized
    file cluster 

    publishDir params.output, saveAs: { filename -> getOutput(params.runid, "pasolli/clusters", filename) }

    output:
    file 'final_clusters.tsv' 
    tuple file(".command.sh"), file(".command.out"), file(".command.err"), file(".command.log")

    shell:
    '''
    cp !{cluster} final_clusters.tsv
    '''
}


def mapJoin(channel_a, channel_b, key_a, key_b){
    channel_a \
        | map{ it -> [it[key_a], it] } \
        | cross(channel_b | map{it -> [it[key_b], it]}) \
        | map { it[0][1] + it[1][1] }
}


/*
* List all files and converts them to tuples.
*/
def collectFiles(dir, sra){
   def fileList = [];
   dir.eachFileRecurse { item ->
        fileList.add([sra, item]);
  }
  return fileList;
}

/**
*
* Dereplicate genomes that are listed in a tsv file containing the columns
* "BIN_ID", "COMPLETENESS", "COVERAGE", "CONTAMINATION", "HETEROGENEITY", "PATH", "N50".
* 
* BIN_ID is a unique identifier for the bin.
* PATH represents either a URL, S3 or local path to a file.
*
**/
workflow wDereplicateFile {
   take:
     genomesTableFile    
   main:
     genomesTableFile | _wDereplicate
   emit:
     _wDereplicate.out
}


/**
*
* This entry point takes the same information as the wDereplicateFile entrypoint,
* with the exception that a channel instead of a tsv file is provided.
*
**/
workflow wDereplicateList {
   take:
     genomes_list 
   main:
    genomes_list | map( it -> "${it['BIN_ID']}\t${it['COMPLETENESS']}\t${it['COVERAGE']}\t${it['CONTAMINATION']}\t${it['HETEROGENEITY']}\t${it['PATH']}\t${it['N50']}" ) \
       | collectFile(seed: "BIN_ID\tCOMPLETENESS\tCOVERAGE\tCONTAMINATION\tHETEROGENEITY\tPATH\tN50", newLine: true, keepHeader: false) \
       | _wDereplicate
   emit:
     _wDereplicate.out
}



workflow _wDereplicate {
   take:
     genomesTableFile
   main:
     defaultMashBuffer = params?.steps?.dereplication?.pasolli?.mashBuffer ?: 500

     genomesTableFile | splitCsv(sep: '\t', header: true) | set { genomesTable }

     // filter genomes by contamination and completeness
     genomesTable | filter({ it.COMPLETENESS.toFloat() >= params?.steps?.dereplication?.pasolli?.minimumCompleteness }) \
       | filter({ it.CONTAMINATION.toFloat() <= params?.steps?.dereplication?.pasolli?.maximumContamination }) \
       | map { line -> [file(line.PATH), line.BIN_ID] } | pMashSketchGenome

     // concatenate (paste) multiple sketches in parallel and compute distance
     pMashSketchGenome.out.sketch | buffer(size: defaultMashBuffer, remainder: true) \
         | pMashPaste | toList() | pMashDist | pClusterDistances | set { clusters }

     pMashSketchGenome.out.stagedGenome | set { stagedGenomeBinIDMapping }

     // select representatives
     pSelectRepresentative(genomesTableFile, clusters) | set { representatives }

     // Check if there are representatives to compare 
     representatives.representatives | splitCsv(sep: '\t') | ifEmpty('DONE') | branch { finalize: it=='DONE' 
         other: it!='DONE' } | set { representativesToCompareC }
     representativesToCompareC.other | multiMap { mags ->
        mag1: mags[0]
        mag2: mags[1]
     } |  set { result }

     // Buffer Mags for ANI comparison and run ANI tool 
     defaultANIBuffer = params?.steps?.dereplication?.pasolli?.ANIBuffer ?: 20
     result.mag1 | map(it -> file(it)) | buffer(size: defaultANIBuffer, remainder: true) | set { mag1 }
     result.mag2 | map(it -> file(it)) | buffer(size: defaultANIBuffer, remainder: true) | set { mag2 }
     pANIb(mag1, mag2)
     pTETRA(mag1, mag2)

     // Prepare output and collect representatives as channel
     pANIb.out | mix(pTETRA.out)  \
       | collectFile(newLine: true)  | splitCsv(sep: '\t', header:false, skip: 0) | set {aniComparisons}

     // Extract BIN_ID for final clustering
     GENOME_A_IDX = 0
     GENOME_B_COMPARISON_IDX = 1
     GENOME_A_ID_IDX = 1
     ANI_RESULT = 2
     stagedGenomeBinIDMapping | cross(aniComparisons) \
        | map { comparison -> [comparison[GENOME_B_COMPARISON_IDX][GENOME_A_ID_IDX], comparison[GENOME_A_IDX][GENOME_A_ID_IDX], comparison[GENOME_B_COMPARISON_IDX][ANI_RESULT]] } \
        | set {aniComparisonsIntermediate}

     GENOME_B_IDX = 0
     GENOME_A_COMPARISON_IDX = 1
     GENOME_B_ID_IDX = 1
     GENOME_A_COMPARISON_PATH = 1
     ANI_RESULT = 2
     stagedGenomeBinIDMapping | cross(aniComparisonsIntermediate) \
        | map { comparison -> [comparison[GENOME_A_COMPARISON_IDX][GENOME_A_COMPARISON_PATH], comparison[GENOME_B_IDX][GENOME_B_ID_IDX], comparison[GENOME_A_COMPARISON_IDX][ANI_RESULT]] } \
        | collectFile(newLine: true){ line -> line.join('\t') } |  set {aniComparisonsFinal}

     pGetCluster(representatives.clusters, aniComparisonsFinal, genomesTableFile)
     pFinalize(representativesToCompareC.finalize, representatives.clusters)
   
     IS_REPRESENTATIVE = 1
     PATH_IDX = 1
     pGetCluster.out.final_clusters | mix(pFinalize.out) | splitCsv(sep: '\t', header: true) \
       | filter({ it.REPRESENTATIVE.toFloat() == IS_REPRESENTATIVE }) | map { it -> it['GENOME'] } \
       | join(genomesTable | map{ bin -> [bin.BIN_ID, bin.PATH] }) | map { bin -> bin[PATH_IDX] } \
       | set{representatives}
//| collectFile(newLine: true) | view() 
  emit:
     representatives
}
