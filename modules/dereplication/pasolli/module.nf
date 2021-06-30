nextflow.enable.dsl=2

params?.steps?.dereplication?.pasolli?.minimumCompleteness = 50
params?.steps?.dereplication?.pasolli?.maximumContamination = 5
params?.steps?.dereplication?.pasolli?.cutoff = 0.05
params?.steps?.dereplication?.pasolli?.pyaniParameters = "-m ANIb"
params?.steps?.dereplication?.pasolli?.representativeAniCutoff = 0.95
params?.steps?.dereplication?.pasolli?.method = "ANI"
params?.steps?.dereplication?.pasolli?.buffer = 10

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
    path bins, stageAs: 'input*.txt'

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
    path sketches

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

    input:
    file('distances.tsv')
    file('mapping.tsv') 

    output:
    file 'distances.mapped.tsv' 

    shell:
    '''
    join -t$'\t'  -2 1 -1 1  <(sort -k 1,1 mapping.tsv)  <(sort -k 1,1 distances.tsv)  | cut -f 2- > distances.col1.mapped.tsv
    join -t$'\t'  -1 1 -2 2   <(sort -k 1,1 mapping.tsv)  <(sort -k 2,2 distances.col1.mapped.tsv) | cut -f 2- > distances.mapped.tsv
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

workflow wDereplicateFile {
   take:
     genomes_table_file    
   main:
     genomes_table_file | _wDereplicate
}


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
     genomes_table_file
   main:
     buffer = params?.steps?.dereplication?.pasolli?.buffer ?: 20

     genomes_table_file | splitCsv(sep: '\t', header: true)  \
       | filter({ it.COMPLETENESS.toFloat() >= params?.steps?.dereplication?.pasolli?.minimumCompleteness }) \
       | filter({ it.CONTAMINATION.toFloat() <= params?.steps?.dereplication?.pasolli?.maximumContamination }) \
       | map { it -> it.PATH } | collect | set {mags} 

     MASH_BUFFER = 5000
     mags | flatten | buffer(size: 5000, remainder: true) | pMashSketch 
     pMashSketch.out.mapping | collectFile(name: 'mapping.txt', newLine: true, skip:0) | set { mappings }
     pMashSketch.out.sketch | collect | pMashDist | set { distances } 
     
     pRenameMashDistances(distances, mappings) | pClusterDistances | set { clusters }

     pSelectRepresentative(genomes_table_file, clusters) | set { representatives }

     representatives.representatives | splitCsv(sep: '\t') | ifEmpty('DONE') | branch { finalize: it=='DONE' 
         other: it!='DONE' } | set { representativesToCompareC }

     representativesToCompareC.other | multiMap { mags ->
        mag1: mags[0]
        mag2: mags[1]
     } |  set { result }

     result.mag1 | map(it -> file(it)) | buffer(size: buffer, remainder: true) | set { mag1 }
     result.mag2 | map(it -> file(it)) | buffer(size: buffer, remainder: true) | set { mag2 }

     pANIb(mag1, mag2)
     pTETRA(mag1, mag2)

     pANIb.out | mix(pTETRA.out) \
       | collectFile(newLine: true) | set { all_ani }

     pGetCluster(representatives.clusters, all_ani, genomes_table_file)
     pFinalize(representativesToCompareC.finalize, representatives.clusters)
   
     pGetCluster.out.final_clusters | splitCsv(sep: '\t', header: true) \
       | filter({ it.REPRESENTATIVE.toFloat() == 1 }) | map { it -> it['GENOME'] } | collectFile(newLine: true) | set { representatives_cluster }

     pFinalize.out | splitCsv(sep: '\t', header: true) \
       | filter({ it.REPRESENTATIVE.toFloat() == 1 }) | map { it -> it['GENOME'] } | collectFile(newLine: true) | set { representatives_finalize }

     representatives_cluster | mix(representatives_finalize) | set{representatives}


  emit:
     representatives
}
