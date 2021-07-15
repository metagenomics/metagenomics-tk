nextflow.enable.dsl=2

process pMash {

    container "quay.io/biocontainers/mash:${params.mash_tag}"

    errorStrategy 'retry'

    label 'large'

    input:
    path bins, stageAs: 'input*.txt'

    output:
    tuple file('distance.tsv'), file('mapping.tsv')

    shell:
    '''
    ls -1 input*.txt > list.txt
    mash sketch -o reference -l list.txt
    mash dist -p !{task.cpus} reference.msh  reference.msh | cut -f 1,2,3 > distance.tsv
    for i in $(ls input*.txt); do echo "$i\t$(readlink -f $i)";  done > mapping.tsv
    '''
}


process pMashSketch {

    container "quay.io/biocontainers/mash:${params.mash_tag}"

    errorStrategy 'retry'

    label 'small'

    input:
    val(num)
    path("input_${num}_*.txt")

    output:
    path('reference.msh'), emit: sketch
    path('mapping.tsv'), emit: mapping

    shell:
    '''
    ls -1 input*.txt > list.txt
    mash sketch -o reference -l list.txt
    for i in $(ls input*.txt); do echo "$i\t$(readlink -f $i)";  done > mapping.tsv
    '''
}


process pMashDist {

    container "quay.io/biocontainers/mash:${params.mash_tag}"

    errorStrategy 'retry'

    label 'large'

    input:
    path sketches, stageAs: 'sketch*.msh'

    output:
    file('distances.tsv')

    shell:
    '''
    mash paste final_sketch !{sketches}
    mash dist -p !{task.cpus} final_sketch.msh final_sketch.msh | cut -f 1,2,3 > distances.tsv
    '''
}

process pRenameMashDistances {

    errorStrategy 'retry'

    label 'tiny'

    input:
    file('distances.tsv')
    file('mapping.tsv') 

    output:
    file 'distances.mapped.tsv' 

    shell:
    '''
    mkdir sortTemp
    join -t$'\t'  -2 1 -1 1  <(sort -T sortTemp -k 1,1 mapping.tsv)  <(sort -T sortTemp -k 1,1 distances.tsv)  | cut -f 2- > distances.col1.mapped.tsv
    join -t$'\t'  -1 1 -2 2   <(sort -T sortTemp -k 1,1 mapping.tsv)  <(sort -T sortTemp -k 2,2 distances.col1.mapped.tsv) | cut -f 2- > distances.mapped.tsv
    '''
}

process pClusterDistances {

    errorStrategy 'retry'

    input:
    file('distances.tsv')

    container "pbelmann/python-env:${params.python_env_tag}"

    label 'medium'

    output:
    tuple file("distances.tsv"), file('out/clusters.tsv') 

    shell:
    '''
    mkdir out
    cluster.py -i distances.tsv -c !{params.steps.dereplication.pasolli.cutoff} -o out
    '''
}

process pSelectRepresentative {

    errorStrategy 'retry'

    input:
    path genome_table
    tuple file("distance"), file("cluster")

    container "pbelmann/python-env:${params.python_env_tag}"

    publishDir "${params.output}/dereplication"

    label 'medium'

    output:
    path("intermediate_clusters.tsv"), emit: clusters
    path("refinement/representatives_to_compare.tsv"), emit: representatives

    shell:
    template 'selectRepresentative.sh'
}


process pANIb {

    errorStrategy 'retry'

    label 'small'

    input:
    file("genome1*") 
    file("genome2*") 

    container "leightonpritchard/average_nucleotide_identity:${params.ani_tag}"

    when:
    params.steps.dereplication.pasolli.method.contains("ANI")

    output:
    file("*.out/out.tsv") 

    shell:
    template 'anib.sh'
}


process pTETRA {

    errorStrategy 'retry'

    label 'small'

    input:
    file("genome1*") 
    file("genome2*") 

    container "leightonpritchard/average_nucleotide_identity:${params.ani_tag}"

    output:
    file("*.out/out.tsv") 

    when:
    params.steps.dereplication.pasolli.method.contains("TETRA")

    shell:
    template 'tetra.sh'
}

process pGetCluster {

    errorStrategy 'retry'

    label 'tiny'

    publishDir "${params.output}/dereplication"

    container "pbelmann/python-env:${params.python_env_tag}"

    input:
    path cluster
    path ani_values
    path genomeAttributes

    output:
    path 'final_clusters.tsv', emit: final_clusters
    path 'ani_values.tsv', emit: ani_values

    shell:
    '''
    mkdir out
    cat <(echo "GENOME_A\tGENOME_B\tANI")  <(sed 's/ /\t/g' !{ani_values}) > ani_values.tsv
    get_final_cluster.py -i !{genomeAttributes} -c !{cluster} -r ani_values.tsv -o out  -a !{params.steps.dereplication.pasolli.representativeAniCutoff}
    cp out/representatives.tsv final_clusters.tsv
    '''
}

process pFinalize {

    errorStrategy 'retry'

    input:
    val finalized
    file cluster 

    publishDir "${params.output}/dereplication"

    output:
    file 'final_clusters.tsv' 

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
* This entry point is highly experimental and should be used for retrieving input data for the dereplication module.
*
**/
workflow wDereplicatePath {
    def baseDir = params.baseDir
    def runID = params.runid

    Channel.from(file(baseDir).list()) | filter({ path -> !(path ==~ /.*summary$/)})  | set { sraDatasets } 
    sraDatasets | map { sra ->  [sra, baseDir + "/" + sra + "/" + runID + "/" ]} \
       | set {sraIDs}

    sraIDs | flatMap { sraID -> collectFiles(file(sraID[1]), sraID[0])} | set {sraFiles}
    sraFiles | filter({ it -> (it[1] ==~ /.*\/binning\/0.1.0\/metabat\/.*.fa$/)}) \
        | map{ sra,f -> [SAMPLE:sra, PATH: baseDir.startsWith("s3://")? "s3:/" + f: f, BIN_ID:file(f).name] } | set{bins}

     sraFiles | filter({ sra, path -> (path ==~ /.*\/magAttributes\/0.2.0\/checkm\/.*.tsv$/)}) \
       | splitCsv(header: ["PATH", "SAMPLE", "BIN_ID", "Marker lineage", "# genomes", "# markers", "# marker sets", "0", "1", "2", "3", "4", "5+", "COMPLETENESS", "CONTAMINATION", "HETEROGENEITY"], sep: '\t') \
       | map { sra, bins -> bins}  \
       | set { checkm }

     sraFiles | filter({ sra, path -> (path ==~ /.*\/binning\/0.1.0\/metabat\/.*_bins_stats.tsv$/)}) \
        | splitCsv(header: true, sep: '\t') | map { sra, bins -> bins} | set{binStats}

     mapJoin(binStats, checkm, "BIN_ID", "BIN_ID") | set {checkmBinStats} 
     mapJoin(checkmBinStats, bins, "file", "BIN_ID") | map( it -> "${it['BIN_ID']}\t${it['COMPLETENESS']}\t${it['COVERAGE']}\t${it['CONTAMINATION']}\t${it['HETEROGENEITY']}\t${it['PATH']}\t${it['N50']}" ) \
         | collectFile(seed: "BIN_ID\tCOMPLETENESS\tCOVERAGE\tCONTAMINATION\tHETEROGENEITY\tPATH\tN50", newLine: true, keepHeader: false) | view() | wDereplicateFile 
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
     // filter table based on provided parameters
     genomesTableFile | splitCsv(sep: '\t', header: true)  \
       | filter({ it.COMPLETENESS.toFloat() >= params?.steps?.dereplication?.pasolli?.minimumCompleteness }) \
       | filter({ it.CONTAMINATION.toFloat() <= params?.steps?.dereplication?.pasolli?.maximumContamination }) \
       | map { it -> it.PATH } | collect | set {mags} 

     //Create Mash buffers 
     defaultMashBuffer = params?.steps?.dereplication?.pasolli?.mashBuffer ?: 500
     mags | flatten | buffer(size: defaultMashBuffer, remainder: true) | set {sketchBuffer}
     sketchBuffer | count() | map{ maxBatchNumber -> [1..maxBatchNumber]} | flatten() | set { sketchBufferCount}

     // Compute mash distance 
     pMashSketch(sketchBufferCount, sketchBuffer)
     pMashSketch.out.mapping | collectFile(name: 'mapping.txt', newLine: true, skip:0) | set { mappings }
     pMashSketch.out.sketch | collect | pMashDist | set { distances } 
     pRenameMashDistances(distances, mappings) | pClusterDistances | set { clusters }

     // Select representatives for every cluster
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
     pANIb.out | mix(pTETRA.out) \
       | collectFile(newLine: true) | set { allAni }

     pGetCluster(representatives.clusters, allAni, genomesTableFile)
     pFinalize(representativesToCompareC.finalize, representatives.clusters)
   
     IS_REPRESENTATIVE = 1
     pGetCluster.out.final_clusters | splitCsv(sep: '\t', header: true) \
       | filter({ it.REPRESENTATIVE.toFloat() == IS_REPRESENTATIVE }) | map { it -> it['GENOME'] } \
       | collectFile(newLine: true) | set { representativesCluster }

     pFinalize.out | splitCsv(sep: '\t', header: true) \
       | filter({ it.REPRESENTATIVE.toFloat() == IS_REPRESENTATIVE }) | map { it -> it['GENOME'] } \
       | collectFile(newLine: true) | set { representativesFinalize }

     representativesCluster | mix(representativesFinalize) | set{representatives}
  emit:
     representatives
}
