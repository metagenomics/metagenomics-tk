nextflow.enable.dsl=2
mode = 'mode_not_set'

def getOutput(SAMPLE, RUNID, TOOL, filename){
    return SAMPLE + '/' + RUNID + '/' + params.modules.annotation.name + '/' + 
          params.modules.annotation.version.major + "." + 
          params.modules.annotation.version.minor + "." + 
          params.modules.annotation.version.patch +
          '/' + TOOL + '/' + filename
}

def set_mode(newMode){
    mode = newMode
}

/**
*
* Diamond is used to search for big input queries in large databases.
* Though not as precise as blast it is way faster if you handle large amounts of data.
* You need to call (and fill out) the aws credential file with -c to use this module! 
*
**/
process pDiamond {
   
      container "${params.diamond_image}"
      
      // Databases will be downloaded to a fixed place so that they can be used by future processes.
      // These fixed place has to be outside of the working-directory to be easy to find for every process.
      // Therefore this place has to be mounted to the docker container to be accessible during run time.
      // Another mount flag is used to get a key file (aws format) into the docker-container. 
      // This file is then used by s5cmd. 
      containerOptions " --user 1000:1000 --volume ${params.databases}:${params.databases} \
      --volume ${params.steps.annotation.s5cmd.keyfile}:/.aws/credentials"
 
      tag "$sample"

      label 'large'

      publishDir params.output, saveAs: { filename -> getOutput("${sample}", params.runid, "diamond", filename) }
      
      // UID mapping does not work for some reason. Every time a database directory is created while running docker,
      // the permissions are set to root. This leads to crashes later on.
      // beforeScript is one way to create a directory outside of Docker to tackle this problem. 
      beforeScript "mkdir -p ${params.databases}"

   input:
      tuple val(sample), file(fasta)
   
   output:
      tuple val("${sample}"), path("diamond.${sample}.out"), emit: results
      tuple file(".command.sh"), file(".command.out"), file(".command.err"), file(".command.log"), emit: log

   shell:
      if( mode == 'S3')
         '''
         # Download the database if there is a more recent one online, or if the size differs.
         # The destination folder should be outside of the working directory to share the database with future processes.
         s5cmd !{params.steps.annotation.s5cmd.params} cp -u -s !{params.steps.annotation.diamond.database} !{params.databases}
         diamond !{params.steps.annotation.diamond.params} --threads !{task.cpus}  --out diamond.!{sample}.out \
         --db !{params.databases}$(basename !{params.steps.annotation.diamond.database}) --query !{fasta}
         '''
      else if( mode == 'local')
         '''
         diamond !{params.steps.annotation.diamond.params} --threads !{task.cpus} --out diamond.!{sample}.out \
         --db !{params.steps.annotation.diamond.database} --query !{fasta}
         '''
      else
         error "Invalid annotation database mode: ${mode}"
}

/**
*
* Prodigal is used to predict genes in fasta files.
*
**/
process pProdigal {

      tag "$sample"

      label 'small'

      container "${params.prodigal_image}"

      publishDir params.output, saveAs: { filename -> getOutput("${sample}", params.runid, "prodigal", filename) }

   input:
      tuple val(sample), file(fasta)

   output:
      tuple val("${sample}"), path("${sample}_prodigal.faa"), emit: amino
      tuple val("${sample}"), path("${sample}_prodigal.fna"), emit: nucleotide
      tuple val("${sample}"), path("${sample}_prodigal.gff"), emit: gff
      tuple file(".command.sh"), file(".command.out"), file(".command.err"), file(".command.log"), emit: log

   shell:
      '''
      prodigal -i !{fasta} -p !{params.steps.annotation.prodigal.mode} -f gff \
      -o !{sample}_prodigal.gff -a !{sample}_prodigal.faa -d !{sample}_prodigal.fna
      '''
}


/**
*
* pKEGGFromDiamond is build to handle results of a diamond/blast search in the outfmt 6 standard.
* These search results will be compared to kegg link-files to produce a file where all available kegg information
* for these search results is collected in a centralized way/file.
* You need to call (and fill out) the aws credential file with -c to use this module!
*
**/
process pKEGGFromDiamond {

      tag "$sample"

      label 'small'

      container "${params.python_env_image}"

      // Databases will be downloaded to a fixed place so that they can be used by future processes.
      // These fixed place has to be outside of the working-directory to be easy to find for every process.
      // Therefore this place has to be mounted to the docker container to be accessible during runtime.
      // Another mount flag is used to get a key file (aws format) into the docker-container. 
      // This file is then used by s5cmd. 
      containerOptions " --user 1000:1000 --volume ${params.databases}:${params.databases} \
      --volume ${params.steps.annotation.s5cmd.keyfile}:/.aws/credentials"

      publishDir params.output, saveAs: { filename -> getOutput("${sample}", params.runid, "keggFromDiamond", filename) }

      // UID mapping does not work for some reason. Every time a database directory is created while running docker,
      // the permissions are set to root. This leads to crashes later on.
      // beforeScript is one way to create a directory outside of Docker to tackle this problem.
      beforeScript "mkdir -p ${params.databases}"

   input:
      tuple val(sample), file(diamond_result)

   output:
      tuple val("${sample}"), path("${sample}_kegg.tsv"), emit: kegg_diamond
      tuple file(".command.sh"), file(".command.out"), file(".command.err"), file(".command.log"), emit: log

   shell:
      if( mode == 'S3')
         '''
         # Download the database if there is a more recent one online, or if the size differs.
         # The destination folder should be outside of the working directory to share the database with future processes.
         s5cmd !{params.steps.annotation.s5cmd.params} cp -u -s !{params.steps.annotation.kegg.database}/* !{params.databases}kegg
         diamond2kegg.py !{diamond_result} !{params.databases}kegg !{sample}_kegg.tsv 
         '''
      else if(mode == 'local')
         '''
         diamond2kegg.py !{diamond_result} !{params.steps.annotation.kegg.database} !{sample}_kegg.tsv
         '''
      else 
         error "Invalid annotation database mode: ${mode}"
}

/**
*
* This entry point uses a file to grab all fasta files in the referenced directories for annotation.
* These fastas will be combined to one file to start the main annotation workflow.
* You need to call (and fill out) the aws credential file with -c to use this module!
* The .tsv file has to have: DATASET PATH entries.
* 
* The "database_mode" is used to choose which database path is expected.
* Use "S3" for object storage based databases, or "local" for offline stored copys.
*
**/
workflow wAnnotateFile {

   take:
      projectTableFile
      database_mode
   main:
      set_mode(database_mode)
      projectTableFile | splitCsv(sep: '\t', header: true) | map{ it -> [it.DATASET, file(it.PATH)] } \
      | collectFile() | map{ it -> [it.name, it]} | _wAnnotation
}


/**
*
* The main annotation workflow. 
* It is build to handle one big input fasta file.
* On this file genes will be predicted, these will be diamond-blasted against kegg. 
* At the end kegg-infos of the results will be collected and presented.
*
**/ 
workflow _wAnnotation {

   take:
      fasta
   main:
      fasta | pProdigal
      pDiamond(pProdigal.out.amino)
      pKEGGFromDiamond(pDiamond.out.results)
      pKEGGFromDiamond.out.kegg_diamond | set { keggAnnotation }
   emit:
      keggAnnotation

   
}
