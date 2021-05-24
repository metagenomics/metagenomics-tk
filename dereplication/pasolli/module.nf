nextflow.enable.dsl=2

params.minimum_completeness = 50
params.maximum_contamination = 5
params.cutoff = 0.05
params.pyani_parameters = "-m ANIb"
params.representative_ani_cutoff = 0.95
params.method = "ANI"
params.buffer = 10

process runMash {

    container "quay.io/biocontainers/mash:${params.mash_tag}"

    errorStrategy 'retry'

    label 'large'

    input:
    path fof, stageAs: 'input*.txt' 

    output:
    tuple file('distance.tsv'), file('mapping.tsv')

    shell:
    '''
    ls -1 input*.txt > input_list.txt
    mash sketch -o reference -l input_list.txt
    mash dist -p !{task.cpus} reference.msh  reference.msh | cut -f 1,2,3 > distance.tsv
    for i in $(ls input*.txt); do echo "$i\t$(readlink -f $i)";  done > mapping.tsv
    '''
}

process renameMashDistances {

    errorStrategy 'retry'

    input:
    tuple file('distances.tsv'), file('mapping.tsv') 

    output:
    file 'distances.mapped.tsv' 

    shell:
    '''
    join -t$'\t'  -2 1 -1 1  <(sort -k 1,1 mapping.tsv)  <(sort -k 1,1 distances.tsv)  | cut -f 2- > distances.col1.mapped.tsv
    join -t$'\t'  -1 1 -2 2   <(sort -k 1,1 mapping.tsv)  <(sort -k 2,2 distances.col1.mapped.tsv) | cut -f 2- > distances.mapped.tsv
    '''
}

process clusterDistances {

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
    cluster.py -i distances.tsv -c !{params.cutoff} -o out
    '''
}

process selectRepresentative {

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


process runANIb {

    errorStrategy 'retry'

    label 'small'

    input:
    file("genome1*") 
    file("genome2*") 

    container "leightonpritchard/average_nucleotide_identity:${params.ani_tag}"

    when:
    params.method.contains("ANI")

    output:
    file("*.out/out.tsv") 

    shell:
    template 'anib.sh'
}


process runTETRA {

    errorStrategy 'retry'

    label 'small'

    input:
    file("genome1*") 
    file("genome2*") 

    container "leightonpritchard/average_nucleotide_identity:${params.ani_tag}"

    output:
    file("*.out/out.tsv") 

    when:
    params.method.contains("TETRA")

    shell:
    template 'tetra.sh'
}

process getCluster {

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
    get_final_cluster.py -i !{genomeAttributes} -c !{cluster} -r ani_values.tsv -o out  -a !{params.dereplication.representative_ani_cutoff}
    cp out/representatives.tsv final_clusters.tsv
    '''
}

process finalize {

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

workflow dereplicate_file {
   take:
     genomes_table_file    
   main:
//     if( !params.output || !params.genomes_table) error "Please provide --genomes_table and --output parameters"
     genomes_table_file | dereplicate
     
}


workflow dereplicate_list {
   take:
     genomes_list    
   main:
     genomes_list | map( it ->  "BIN_ID\tCOMPLETENESS\tCOVERAGE\tCONTAMINATION\tHETEROGENEITY\tPATH\tN50\n${it['BIN_ID']}\t${it['COMPLETENESS']}\t${it['COVERAGE']}\t${it['CONTAMINATION']}\t${it['HETEROGENEITY']}\t${it['PATH']}\t${it['N50']}\n" ) | collectFile(newLine: false, keepHeader: true) | dereplicate 
   emit:
     dereplicate.out
}


workflow dereplicate {
   take:
     genomes_table_file
   main:
     genomes_table_file | splitCsv(sep: '\t', header: true)  \
       | filter({ it.COMPLETENESS.toFloat() >= params.dereplication.minimum_completeness }) \
       | filter({ it.CONTAMINATION.toFloat() <= params.dereplication.maximum_contamination }) \
       | map { it -> it.PATH } | collect | set {mags} 

     mags | runMash | renameMashDistances | clusterDistances | set { clusters }

     selectRepresentative(genomes_table_file, clusters) | set { representatives }

     representatives.representatives | splitCsv(sep: '\t') | ifEmpty('DONE') | branch { finalize: it=='DONE' 
         other: it!='DONE' } | set { representativesToCompareC }

     representativesToCompareC.other | multiMap { mags ->
        mag1: mags[0]
        mag2: mags[1]
     } |  set { result }

     result.mag1 | map(it -> file(it)) | buffer(size: params.buffer, remainder: true) | set { mag1 }
     result.mag2 | map(it -> file(it)) | buffer(size: params.buffer, remainder: true) | set { mag2 }

     runANIb(mag1, mag2)
     runTETRA(mag1, mag2)

     runANIb.out | mix(runTETRA.out) \
       | collectFile(newLine: true) | set { all_ani }

     getCluster(representatives.clusters, all_ani, genomes_table_file)
     finalize(representativesToCompareC.finalize, representatives.clusters)
   
     getCluster.out.final_clusters | splitCsv(sep: '\t', header: true) \
       | filter({ it.REPRESENTATIVE.toFloat() == 1 }) | map { it -> it['GENOME'] } | collectFile(newLine: true) | set { representatives_cluster }

     finalize.out | splitCsv(sep: '\t', header: true) \
       | filter({ it.REPRESENTATIVE.toFloat() == 1 }) | map { it -> it['GENOME'] } | collectFile(newLine: true) | set { representatives_finalize }

     representatives_cluster | mix(representatives_finalize) | set{representatives}


  emit:
     representatives
}


workflow {
   dereplicate_file(Channel.value(file(params.genomes_table)))
}
