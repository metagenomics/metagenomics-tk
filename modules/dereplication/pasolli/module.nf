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
    mash sketch !{binid}  -o !{binid}.msh
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
    mash dist -p !{task.cpus} final_sketch.msh final_sketch.msh | cut -f 1,2,3 > distances.tsv
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
    cluster.py -i distances.tsv -c !{params.steps.dereplication.pasolli.cutoff} -o out
    '''
}


process pSelectRepresentative {

    errorStrategy 'retry'

    input:
    path genome_table
    tuple file("distance"), file("cluster")

    container "${params.python_env_image}"

    publishDir params.output, saveAs: { filename -> getOutput(params.runid, "pasolli/species/selectedRepresentatives", filename) }

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

    publishDir params.output, saveAs: { filename -> getOutput(params.runid, "pasolli/species/clusters", filename) }

    container "${params.python_env_image}"

    input:
    path cluster, stageAs: 'cluster'
    path ani_values, stageAs: 'ani_values'
    path genomeAttributes, stageAs: 'genomeAttributes'

    output:
    path 'clusters.tsv', emit: finalClusters
    path 'ani_values.tsv', emit: ani_values
    tuple file(".command.sh"), file(".command.out"), file(".command.err"), file(".command.log")

    shell:
    '''
    mkdir out
    cat <(echo "GENOME_A\tGENOME_B\tANI")  <(sed 's/ /\t/g' !{ani_values}) > ani_values.tsv
    get_final_cluster.py -i !{genomeAttributes} -c !{cluster} -r ani_values.tsv -o . -a !{params.steps.dereplication.pasolli.representativeAniCutoff}
    '''
}


process pFinalize {

    input:
    val finalized
    file cluster 

    publishDir params.output, saveAs: { filename -> getOutput(params.runid, "pasolli/species/clusters", filename) }

    output:
    file 'clusters.tsv' 
    tuple file(".command.sh"), file(".command.out"), file(".command.err"), file(".command.log")

    shell:
    '''
    cp !{cluster} clusters.tsv
    '''
}


process pSANS {

    container "${params.sans_image}"

    label 'small'

    tag "Cluster ${clusterID}"

    publishDir params.output, saveAs: { filename -> getOutput(params.runid, "pasolli/strains/sans", filename) }

    when:
    params.steps.dereplication.containsKey("sans")

    input:
    path(ids)
    path(genomes, stageAs: 'genome_?') 
    val(clusterID)

    output:
    tuple val("${clusterID}"), file("${clusterID}_newick.txt"), emit: newick
    tuple val("${clusterID}"), file("${clusterID}_clusters.tsv"), emit: clusters

    shell:
    '''
    mkdir input
    # get ID from list of genomes
    for g in genome_* ; do 
	line=$(echo $g | cut -d '_' -f 2)
	ID=$(sed "${line}q;d" !{ids})
	ln -s $(readlink -f $g) input/${ID}
    done 
    # Run SANS for a set of genomes and translate the newick tree to a set of clusters
    find $PWD/input -type l > input.tsv
    /sans/SANS-autoN.sh -i input.tsv !{params.steps.dereplication.sans.additionalParams} -N !{clusterID}_newick.txt
    /sans/scripts/newick2clusters.py !{clusterID}_newick.txt > !{clusterID}_clusters.tsv
    '''
}


process pGetStrainClusterRepresentatives {

    label 'tiny'

    publishDir params.output, saveAs: { filename -> getOutput(params.runid, "pasolli/strains/", filename) }

    container "${params.python_env_image}"

    when:
    params.steps.dereplication.containsKey("sans")

    input:
    path cluster
    path genomeAttributes

    output:
    path 'clusters.tsv', emit: strainRepresentatives

    shell:
    '''
    select_strain_representative.py -i !{genomeAttributes} -c !{cluster} -o .
    '''
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


/**
*
* This workflow accepts species clusters (speciesClusters) and a table (genomesTableFile) containing bin attributes as channel inputs. 
* The species clusters channel has the following format: [BIN_ID, PATH, CLUSTER].
* The genome table file contains columns such as N50, COMPLETENESS and CONTAMINATION.
*/
workflow _wStrainDereplication {
   take:
     speciesClusters
     genomesTableFile
   main:
     CLUSTER_NAME_IDX=2
     BIN_ID_IDX=0
    
     // Species clusters consisting of just one genome do not have to be dereplicated again.
     speciesClusters | groupTuple(by: CLUSTER_NAME_IDX) | set { groupedSpecies } 
     groupedSpecies | collectFile(newLine: true){ ids -> [ids[CLUSTER_NAME_IDX], ids[BIN_ID_IDX].join('\n')] } \
         | map { cluster -> [cluster.name, cluster] } | set { binIDs }

     // Species clusters consisting of just of genome do not have to be dereplicated again.
     BINS_IDX=2
     CLUSTER_INDEX = 0
     groupedSpecies | map { it -> it.reverse(false) }  | join(binIDs, by: CLUSTER_INDEX) \
      | groupTuple(by: CLUSTER_NAME_IDX) \
      | branch { 
        single: it[BINS_IDX].size() == 1
        multi: it[BINS_IDX].size() > 1
     } | set { clusterSize }  

     // Split values of an array (BIN ID, BIN PATH and CLUSTER ID) in different channels, that can
     // be consumed by the SANS process
     ARRAY_CLUSTER_IDX=0
     ARRAY_BIN_PATH_IDX=1
     ARRAY_BIN_ID_IDX=3
     clusterSize.multi  | multiMap { mags ->
        binIds: mags[ARRAY_BIN_ID_IDX]
        binPaths: mags[ARRAY_BIN_PATH_IDX][0]
        clusterId: mags[ARRAY_CLUSTER_IDX].flatten()
     } | set { result }

     pSANS(result.binIds, result.binPaths , result.clusterId | flatten)

     // Combine all SANS cluster and select representatives
     STRAIN_IDX=1
     STRAIN_CLUSTER_IDX=1
     STRAIN_BIN_PATH_IDX=0
     SPECIES_CLUSTER_IDX=0
     pSANS.out.clusters | splitCsv(sep: '\t') \
	| map { it -> [it[SPECIES_CLUSTER_IDX], file(it[STRAIN_IDX][STRAIN_BIN_PATH_IDX]).name,it[STRAIN_IDX][STRAIN_CLUSTER_IDX]] } \
	| set {SANSCluster}

     CLUSTER_IDX=0
     BIN_ID_LABEL_IDX=2
     SINGLE_CLUSTER_ID=0
     clusterSize.single | map { it -> [it[CLUSTER_IDX][CLUSTER_IDX], *it[BIN_ID_LABEL_IDX], SINGLE_CLUSTER_ID] } \
	| mix(SANSCluster) \
	| collectFile(newLine:true, seed: "CLUSTER\tBIN_ID\tSTRAIN_CLUSTER"){ it -> ["strainCluster.tsv", it.join('\t')] } \
	| set { strainRepresentatives }
     
     pGetStrainClusterRepresentatives(strainRepresentatives, genomesTableFile)
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
       | collectFile(newLine: true) | splitCsv(sep: '\t', header:false, skip: 0) | set {aniComparisons}

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
   
     // Prepare genome files for strain dereplication
     IS_REPRESENTATIVE = 1
     PATH_IDX = 1
     pGetCluster.out.finalClusters | mix(pFinalize.out) | splitCsv(sep: '\t', header: true) \
	| set { finalClusters  }

     finalClusters | map { bin -> [bin.CLUSTER, bin.GENOME] } | set{ clustersGenome } 
     genomesTable | map { bin -> [bin.PATH, bin.BIN_ID] } \
	| join(clustersGenome, by: 1) |  set { clusterFiles  }

     _wStrainDereplication(clusterFiles, genomesTableFile)

     finalClusters | filter({ it.REPRESENTATIVE.toFloat() == IS_REPRESENTATIVE }) | map { it -> it['GENOME'] } \
       | join(genomesTable | map{ bin -> [bin.BIN_ID, bin.PATH] }) | map { bin -> bin[PATH_IDX] } \
       | set{representatives}

  emit:
     representatives
}
