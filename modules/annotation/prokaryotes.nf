include { _wAnalyseMMseqs2; _wAnalyseMMseqs2Taxonomy } from './mmseqs2'
include { _wAnalyseDiamond } from './diamond'
include { pDumpLogs } from '../utils/processes'
include { checkAnnotationSearchTool } from '../utils/methods'

/**
*
* pKEGGFromBlast is build to handle results in the outfmt 6 file standard.
* The results can be produced by MMseqs2 or DIAMOND, only the columns query and target are used.
* These search results will be compared to kegg link-files to produce a file where all available kegg information
* for these search results is collected in a centralized way/file.
* You need to call (and fill out) the aws credential file with -c to use this module as kegg is commercial an the database has to be placed in your project s3!
*
**/
process pKEGGFromBlast {

      tag "$sample"

      label 'small'

      container "${params.python_env_image}"

      // Databases will be downloaded to a fixed place so that they can be used by future processes.
      // These fixed place has to be outside of the working-directory to be easy to find for every process.
      // Therefore this place has to be mounted to the docker container to be accessible during runtime.
      // Another mount flag is used to get a key file (aws format) into the docker-container. 
      // This file is then used by s5cmd. 

      containerOptions Utils.getDockerMount(params.steps?.annotation?.keggFromBlast?.database, params, apptainer=params.apptainer) + (params.apptainer ? "" : Utils.getDockerNetwork())

      publishDir params.output, mode: "${params.publishDirMode}", saveAs: { filename -> Output.getOutput("${sample}", params.runid, "keggFromBlast", params.modules.annotation, filename) }, \
         pattern: "{**.tsv}"

      when params?.steps.containsKey("annotation") && params?.steps.annotation.containsKey("keggFromBlast")

   secret { "${S3_KEGG_ACCESS}"!="" ? ["S3_kegg_ACCESS", "S3_kegg_SECRET"] : [] } 

   input:
      tuple val(sample), val(binType), file(blastResult)

   output:
      tuple val("${sample}"), path("*.keggPaths.tsv"), emit: keggPaths
      tuple val("${sample}_${binType}"), val("${output}"), val(params.LOG_LEVELS.INFO), file(".command.sh"), \
        file(".command.out"), file(".command.err"), file(".command.log"), emit: logs

   shell:
      output = Output.getOutput("${sample}", params.runid, "keggFromBlast", params.modules.annotation, "")
      DOWNLOAD_LINK=params.steps?.annotation?.keggFromBlast?.database?.download?.source ?: ""
      MD5SUM=params?.steps?.annotation?.keggFromBlast?.database?.download?.md5sum ?: ""
      S5CMD_PARAMS=params.steps?.annotation?.keggFromBlast?.database?.download?.s5cmd?.params ?: ""
      EXTRACTED_DB=params.steps?.annotation?.keggFromBlast?.database?.extractedDBPath ?: ""
      S3_KEGG_ACCESS=params.steps?.annotation?.keggFromBlast?.database?.download?.s5cmd && S5CMD_PARAMS.indexOf("--no-sign-request") == -1 ? "\$S3_kegg_ACCESS" : ""
      S3_KEGG_SECRET=params.steps?.annotation?.keggFromBlast?.database?.download?.s5cmd && S5CMD_PARAMS.indexOf("--no-sign-request") == -1 ? "\$S3_kegg_SECRET" : ""
      '''
      # Check developer documentation
      KEGG_DB=""
      if [[ -z "!{EXTRACTED_DB}" ]] 
      then
        DATABASE=!{params.polished.databases}/kegg
        LOCK_FILE=${DATABASE}/lock.txt

	# Check if access and secret keys are necessary for s5cmd
        if [ ! -z "!{S3_KEGG_ACCESS}" ]
        then
          export AWS_ACCESS_KEY_ID=!{S3_KEGG_ACCESS}
          export AWS_SECRET_ACCESS_KEY=!{S3_KEGG_SECRET}
        fi

        # Download KEGG database
        mkdir -p ${DATABASE}
        flock ${LOCK_FILE} concurrentDownload.sh --output=${DATABASE} \
         --link=!{DOWNLOAD_LINK} \
         --httpsCommand="wgetStatic --no-check-certificate -qO- !{DOWNLOAD_LINK} | tar -xz " \
         --s3DirectoryCommand="s5cmd !{S5CMD_PARAMS} cp --concurrency !{task.cpus} !{DOWNLOAD_LINK} . " \
         --s3FileCommand="s5cmd !{S5CMD_PARAMS} cat !{DOWNLOAD_LINK} | tar -xz " \
	     --s5cmdAdditionalParams="!{S5CMD_PARAMS}" \
         --localCommand="tar -xzvf !{DOWNLOAD_LINK} " \
         --expectedMD5SUM=!{MD5SUM}

         KEGG_DB="${DATABASE}/out/"
      else
         KEGG_DB="!{EXTRACTED_DB}"
      fi
      blast2kegg.py !{blastResult} ${KEGG_DB} !{blastResult.baseName}.keggPaths.tsv
      '''
}


/**
*
* Analyses the predicted proteins of prokaryotic sequences.
* DIAMOND or MMseqs2 is used for protein search.
* The taxonomic classification is done via MMseqs2 (mmseqs2_taxonomy).
*
**/
workflow _wAnalyseProkaryotes {
  take:
    sourceChannel
    contig2GeneMapping
    proteins
    fastaCounter
  main:
    checkAnnotationSearchTool(params?.steps?.annotation ?: [:])

    // The proteins are searched with either DIAMOND or MMseqs2, never with both
    if (params?.steps?.annotation?.containsKey("diamond")) {
      _wAnalyseDiamond(sourceChannel, contig2GeneMapping, proteins, fastaCounter)
      _wAnalyseDiamond.out.blast | set { collectedBlast }
    } else {
      _wAnalyseMMseqs2(sourceChannel, contig2GeneMapping, proteins, fastaCounter)
      _wAnalyseMMseqs2.out.blast | set { collectedBlast }
    }

    // The taxonomic classification is independent of the search tool
    _wAnalyseMMseqs2Taxonomy(sourceChannel, proteins, fastaCounter)

    // The kegg database result is selected and adjusted for further analysis via pKEGGFromBlast
    collectedBlast
        | filter { sample, type, dbType, blastFiles -> dbType == "kegg" }
        | map { sample, type, dbType, blastFiles -> [sample, type, blastFiles] }
        | pKEGGFromBlast

    pKEGGFromBlast.out.logs | pDumpLogs
  emit:
    keggAnnotation = pKEGGFromBlast.out.keggPaths
    mmseqs2_kronaHtml = _wAnalyseMMseqs2Taxonomy.out.mmseqs2_kronaHtml
    mmseqs2_krakenTaxonomy = _wAnalyseMMseqs2Taxonomy.out.mmseqs2_krakenTaxonomy
    mmseqs2_taxonomy = _wAnalyseMMseqs2Taxonomy.out.mmseqs2_taxonomy
    blast = collectedBlast
}
