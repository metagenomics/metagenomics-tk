nextflow.enable.dsl=2

include { pDumpLogs } from '../utils/processes'

mode = 'mode_not_set'

def getOutput(SAMPLE, RUNID, TOOL, filename){
    return SAMPLE + '/' + RUNID + '/' + params.modules.annotation.name + '/' + 
          params.modules.annotation.version.major + "." + 
          params.modules.annotation.version.minor + "." + 
          params.modules.annotation.version.patch +
          '/' + TOOL + '/' + filename
}

/**
* The "database_mode" is used to choose which database path is expected.
* If the Diamond-databasepath starts with "https://" or "s3://" the object storage based mode is used.
* All other paths are seen as "local" mode paths and offline stored copys are expected.
**/
def set_mode(pathString){
    if(pathString != null){
      if (pathString.startsWith("s3://")){
        mode = "S3";
      } else if (pathString.startsWith("https://")) {
        mode = "S3";
      } else {
        mode = "local";
      }
    }
}


>>>>>>> master
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

      publishDir params.output, saveAs: { filename -> getOutput("${sample}", params.runid, "diamond", filename) }, \
         pattern: "{**.diamond.out}"
      
      // UID mapping does not work for some reason. Every time a database directory is created while running docker,
      // the permissions are set to root. This leads to crashes later on.
      // beforeScript is one way to create a directory outside of Docker to tackle this problem. 
      beforeScript "mkdir -p ${params.databases}"

      when params?.steps.containsKey("annotation") && params?.steps.annotation.containsKey("diamond")

   input:
      tuple val(sample), val(binID), file(fasta)
   
   output:
      tuple val("${sample}"), val("${binID}"), path("${binID}.diamond.out"), emit: results
      tuple val("${sample}_${binID}"), val("${output}"), val(params.LOG_LEVELS.INFO), file(".command.sh"), \
        file(".command.out"), file(".command.err"), file(".command.log"), emit: logs


   shell:
      output = getOutput("${sample}", params.runid, "diamond", "")
      if( mode == 'S3')
         '''
         # Download the database if there is a more recent one online, or if the size differs.
         # The destination folder should be outside of the working directory to share the database with future processes.
         DATABASE=!{params.databases}/diamond
         mkdir -p ${DATABASE}
         DIAMOND_PATH=$(readlink -f ${DATABASE}/*)
         s5cmd !{params.steps.annotation.s5cmd.params} cp -u -s !{params.steps.annotation.diamond.database} ${DATABASE}
         diamond !{params.steps.annotation.diamond.params} --threads !{task.cpus}  --out !{binID}.diamond.out \
         --db ${DIAMOND_PATH} --query !{fasta}
         '''
      else if( mode == 'local')
         '''
         diamond !{params.steps.annotation.diamond.params} --threads !{task.cpus} --out !{binID}.diamond.out \
         --db !{params.steps.annotation.diamond.database} --query !{fasta}
         '''
      else
         error "Invalid annotation database mode: ${mode}"
}


process pResistanceGeneIdentifier {
   
      container "${params.rgi_image}"
      
      containerOptions " --user 1000:1000 --volume ${params.databases}:${params.databases} \
                        --volume ${params.steps.annotation.s5cmd.keyfile}:/.aws/credentials"
 
      tag "$sample $binID"

      label 'large'

      publishDir params.output, saveAs: { filename -> getOutput("${sample}", params.runid, "rgi", filename) }, \
         pattern: "{**.rgi.txt}"

      
      beforeScript "mkdir -p ${params.databases}"

      when params.steps.containsKey("annotation") && params?.steps.annotation.containsKey("rgi")

   input:
      tuple val(sample), val(binID), file(fasta)
   
   output:
      tuple val("${sample}"), val("${binID}"), path("${binID}.rgi.txt"), emit: results
      tuple val("${sample}_${binID}"), val("${output}"), val(params.LOG_LEVELS.INFO), file(".command.sh"), \
        file(".command.out"), file(".command.err"), file(".command.log"), emit: logs

   shell:
      output = getOutput("${sample}", params.runid, "rgi", "")
      '''
      DATABASE=!{params.databases}/rgi
      LOCK_FILE=${DATABASE}/checksum.txt
      DOWNLOAD_LINK=!{params.steps.annotation.rgi.source}
      MD5SUM=!{params.steps.annotation.rgi.md5sum}

      mkdir -p ${DATABASE}
      flock ${LOCK_FILE} concurrentDownload.sh --output=${DATABASE} \
         --link=$DOWNLOAD_LINK \
         --httpsCommand="wget -O data ${DOWNLOAD_LINK} && tar -xvf data ./card.json && rm data" \
         --s3Command="s5cmd !{params.steps?.annotation?.s5cmd.params} cp ${DOWNLOAD_LINK} data && tar -xvf data ./card.json && rm data" \
         --localCommand="tar -xvf ${DOWNLOAD_LINK} ./card.json" \
         --expectedMD5SUM=${MD5SUM}

       sed 's/*//g' !{fasta} > input.faa
       rgi load --card_json ${DATABASE}/out/card.json --local
       rgi main --input_sequence input.faa \
                  --output_file !{binID}.rgi --input_type protein --local \
                  --alignment_tool DIAMOND --num_threads !{task.cpus} --clean
       '''
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

      publishDir params.output, saveAs: { filename -> getOutput("${sample}", params.runid, "prodigal", filename) }, \
         pattern: "{**.faa,**.fna,**.gff}"

      when params.steps.containsKey("annotation") && params.steps.annotation.containsKey("prodigal")

   input:
      val(mode)
      tuple val(sample), val(binID), file(fasta)

   output:
      tuple val("${sample}"), val("${binID}"), path("${binID}_prodigal.faa"), emit: amino
      tuple val("${sample}"), val("${binID}"), path("${binID}_prodigal.fna"), emit: nucleotide
      tuple val("${sample}"), val("${binID}"), path("${binID}_prodigal.gff"), emit: gff
      tuple val("${sample}_${binID}"), val("${output}"), val(params.LOG_LEVELS.INFO), file(".command.sh"), \
        file(".command.out"), file(".command.err"), file(".command.log"), emit: logs


   shell:
      output = getOutput("${sample}", params.runid, "prodigal", "")
      if(mode == "param")
        '''
        pigz -f -dc !{fasta} | prodigal -p !{params.steps.annotation.prodigal.mode} -f gff \
        -o !{binID}_prodigal.gff -a !{binID}_prodigal.faa -d !{binID}_prodigal.fna
        '''
      else if(mode == "single")
        '''
        pigz -f -dc !{fasta} | prodigal -p single -f gff \
        -o !{binID}_prodigal.gff -a !{binID}_prodigal.faa -d !{binID}_prodigal.fna
        '''
      else if(mode == "meta")
        '''
        pigz -f -dc !{fasta} | prodigal -p meta -f gff \
        -o !{binID}_prodigal.gff -a !{binID}_prodigal.faa -d !{binID}_prodigal.fna
        '''
      else
         error "Invalid prodigal mode: ${mode}"
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

      publishDir params.output, saveAs: { filename -> getOutput("${sample}", params.runid, "keggFromDiamond", filename) }, \
         pattern: "{**.tsv}"

      // UID mapping does not work for some reason. Every time a database directory is created while running docker,
      // the permissions are set to root. This leads to crashes later on.
      // beforeScript is one way to create a directory outside of Docker to tackle this problem.
      beforeScript "mkdir -p ${params.databases}"
      when params?.steps.containsKey("annotation") && params?.steps.annotation.containsKey("keggFromDiamond")

   input:
      tuple val(sample), val(binID), file(diamond_result)

   output:
      tuple val("${sample}"), path("${binID}_kegg.tsv"), emit: kegg_diamond
      tuple val("${sample}_${binID}"), val("${output}"), val(params.LOG_LEVELS.INFO), file(".command.sh"), \
        file(".command.out"), file(".command.err"), file(".command.log"), emit: logs

   shell:
      output = getOutput("${sample}", params.runid, "keggFromDiamond", "")
      if( mode == 'S3')
         '''
         # Download the database if there is a more recent one online, or if the size differs.
         # The destination folder should be outside of the working directory to share the database with future processes.
         s5cmd !{params.steps.annotation.s5cmd.params} cp -u -s !{params.steps.annotation.kegg.database}/* !{params.databases}kegg
         diamond2kegg.py !{diamond_result} !{params.databases}kegg !{binID}_kegg.tsv 
         '''
      else if(mode == 'local')
         '''
         diamond2kegg.py !{diamond_result} !{params.steps.annotation.kegg.database} !{binID}_kegg.tsv
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
* If the Diamond-databasepath starts with "https://" or "s3://" the object storage based mode is used.
* All other paths are seen as "local" mode paths and offline stored copys are expected.
*
**/
workflow wAnnotateFile {

   take:
      projectTableFile
   main:
      annotationTmpDir = params.tempdir + "/annotation"
      file(annotationTmpDir).mkdirs()
      set_mode(params.steps.annotation.diamond.database)
      projectTableFile | splitCsv(sep: '\t', header: true) \
      | map{ [it.DATASET, file(it.PATH)] } \
      | collectFile(tempDir: params.tempdir + "/annotation") | map{ [it.name, it]} | set { input } 
      | _wAnnotation(Channel.value("param"), input)
   emit:
      keggAnnotation = _wAnnotation.out.keggAnnotation
}


def flattenBins(binning){
  def chunkList = [];
  def SAMPLE_IDX = 0;
  def BIN_PATHS_IDX = 1;
  binning[BIN_PATHS_IDX].each {
     chunkList.add([binning[SAMPLE_IDX], it]);
  }
  return chunkList;
}


/**
*
* See wAnnotateFile for a description.
* This entry point concatenates all files of a sample and annotates these. 
*
**/
workflow wAnnotateList {
   take:
      prodigalMode
      fasta
   main:
      annotationTmpDir = params.tempdir + "/annotation"
      file(annotationTmpDir).mkdirs()
      set_mode(params.steps?.annotation?.diamond?.database)
      _wAnnotation(prodigalMode, fasta)
    emit:
      keggAnnotation = _wAnnotation.out.keggAnnotation
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
      prodigalMode
      fasta
   main:
      pProdigal(prodigalMode,fasta)
      pDiamond(pProdigal.out.amino)
      pProdigal.out.amino | pResistanceGeneIdentifier
      pKEGGFromDiamond(pDiamond.out.results)
      pKEGGFromDiamond.out.kegg_diamond | set { keggAnnotation }

      pProdigal.out.logs | mix(pDiamond.out.logs) \
	| mix(pResistanceGeneIdentifier.out.logs) \
	| mix(pKEGGFromDiamond.out.logs) | pDumpLogs

   emit:
      keggAnnotation   
}
