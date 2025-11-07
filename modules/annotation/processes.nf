/**
*
* The pCollectFile process is designed to collect and concatenate BLAST and other MMseqs/MetaEuck output files.
* It first concatenates all files and splits them based on the BIN_ID.
*
**/
process pCollectFile {

    memory { Utils.getMemoryResources(params.resources.small, "${sample}", task.attempt, params.resources) }

    cpus { Utils.getCPUsResources(params.resources.small, "${sample}", task.attempt, params.resources) }

    tag "Sample: ${sample}, Database: ${dbType}"

    container "${params.ubuntu_image}"

    publishDir params.output, mode: "${params.publishDirMode}", saveAs: { filename -> Output.getOutput("${sample}", params.runid, "${tool}/${dbType}", params.modules.annotation, filename) }

    input:
    tuple val(dbType), val(sample), val(type), path(blastOutputs), val(tool), val(fileType)

    output:
    tuple val("${sample}"), val("${type}"), val("${dbType}"), path("*${sample}_${type}.${dbType}.blast.tsv", arity: "0..*"), optional: true, emit: blast
    tuple val("${sample}"), val("${type}"), val("${dbType}"), path("*${sample}_${type}.${dbType}.tax_per_pred.tsv", arity: "0..*"), optional: true, emit: taxPred
    tuple val("${sample}"), val("${type}"), val("${dbType}"), path("*${sample}_${type}.${dbType}.tax_per_contig.tsv", arity: "0..*"), optional: true, emit: taxPerContig
    tuple val("${sample}"), val("${type}"), val("${dbType}"), path("*${sample}_${type}.${dbType}.gff", arity: "0..*"), optional: true, emit: gff
    tuple val("${sample}"), val("${type}"), val("${dbType}"), path("*${sample}_${type}.${dbType}.faa", arity: "0..*"), optional: true, emit: faa
    tuple val("${sample}"), val("${type}"), val("${dbType}"), path("*${sample}_${type}.${dbType}.fas", arity: "0..*"), optional: true, emit: fas

    shell:
    if (fileType == "blast") {
        '''
    mkdir tmp

    # Concatenate all chunks
    mlr --itsv --otsv cat  *.tsv  > tmp/concat.tsv

    # Create for every bin a seperate file in parallel  
    csvtk cut -t -T  -f BIN_ID tmp/concat.tsv \
	| tail -n +2 | sort | uniq \
	| xargs -P !{task.cpus} -I {} bash filterBin.sh {} "!{sample}_!{type}.!{dbType}.blast.tsv" 
    '''
    }
    else if (fileType == "gff") {
        '''
    mkdir tmp

    # Concatenate all chunks
    mlr --itsv --otsv cat  *.gff  > tmp/concat.gff

    # Create for every bin a seperate file in parallel  
    csvtk cut -t -T  -f BIN_ID tmp/concat.gff \
	| tail -n +2 | sort | uniq \
	| xargs -P !{task.cpus} -I {} bash filterBin.sh {} "!{sample}_!{type}.!{dbType}.gff"
    '''
    }
    else if (fileType == "taxPred") {
        '''
    mkdir tmp

    # Concatenate all chunks
    mlr --itsv --otsv cat  *.tax_per_pred.tsv  > tmp/all.tax_per_pred.tsv

    # Create for every bin a seperate file in parallel  
    csvtk cut -t -T  -f BIN_ID tmp/all.tax_per_pred.tsv \
	| tail -n +2 | sort | uniq \
	| xargs -P !{task.cpus} -I {} bash filterBin.sh {} "!{sample}_!{type}.!{dbType}.tax_per_pred.tsv"
    '''
    }
    else if (fileType == "taxPerContig") {
        '''
    mkdir tmp

    # Concatenate all chunks
    mlr --itsv --otsv cat  *.tax_per_contig.tsv  > tmp/all.tax_per_contig.tsv

    # Create for every bin a seperate file in parallel  
    csvtk cut -t -T  -f BIN_ID tmp/all.tax_per_contig.tsv \
	| tail -n +2 | sort | uniq \
	| xargs -P !{task.cpus} -I {} bash filterBin.sh {} "!{sample}_!{type}.!{dbType}.tax_per_contig.tsv"
    '''
    }
    else if (fileType == "codon") {
        template("filterCodon.sh")
    }
    else if (fileType == "faa") {
        template("filterFaa.sh")
    }
}


workflow _wCollectMetaEuk {
    take:
    metaEukOut
    fileEnding

    main:
    FIRST_ELEM_IDX = 0
    FIRST_ELEM_DB_IDX = 0
    FIRST_ELEM_SAMPLE_IDX = 1
    FIRST_ELEM_TYPE_IDX = 2

    ELEM_PATH_IDX = 3

    metaEukOut
        | map { db, sample, type, start, stop, chunks, out -> [db + "_-_" + sample + "_-_" + type, db, sample, type, start, stop, chunks, out] }
        | map { key, db, sample, type, _start, _stop, chunks, out -> tuple(groupKey(key, chunks.toInteger()), [db, sample, type, out]) }
        | groupTuple()
        | map { _key, dataset ->
            [
                dataset[FIRST_ELEM_IDX][FIRST_ELEM_DB_IDX],
                dataset[FIRST_ELEM_IDX][FIRST_ELEM_SAMPLE_IDX],
                dataset[FIRST_ELEM_IDX][FIRST_ELEM_TYPE_IDX],
                dataset.stream().map { elem -> elem[ELEM_PATH_IDX] }.collect(),
            ]
        }
        | combine(Channel.value("metaeuk"))
        | combine(fileEnding)
        | pCollectFile
}
