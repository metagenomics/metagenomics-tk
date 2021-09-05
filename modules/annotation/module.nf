nextflow.enable.dsl=2

MODULE="annotation"
VERSION="0.1.0"

def getOutput(SAMPLE, RUNID, TOOL, filename){
    return SAMPLE + '/' + RUNID + '/' + MODULE + '/' + VERSION + '/' + TOOL + '/' + filename
}


/**
process pBlast {

   container "ncbi/blast:${params.blast_tag}"

   tag "$sample"

   label 'medium'

   publishDir params.output, saveAs: { filename -> getOutput("${sample}", params.runid, "blast", filename) }

   input:
   //Überlegen, wie man sämtliche S3 Inputs herunterläd und danach parallelisiert auf die Blast Prozesse verteilt.
  // Vieleicht doch zuerst eine Art staging Prozess?

} **/


process pDiamond {

   container "quay.io/biocontainers/diamond:${params.diamond_tag}"

   tag "$sample"

   label 'large'

   publishDir params.output, saveAs: { filename -> getOutput("${sample}", params.runid, "diamond", filename) }

   beforeScript 'echo Test' 

   input:
   tuple val(sample), val(project)
   output:
   tuple val("${sample}"), path("*.daa"), emit: diamond_results

   shell:
   '''
   #s5cmd !{params.steps.annotation.diamond.s5cmd} cp -u !{params.steps.annotation.diamond.database} databases
   mkdir input_diamond
   #s5cmd !{params.steps.annotation.diamond.s5cmd} cp --flatten !{project}/*/binning/*/metabat/*.fa ./input_diamond
   
   echo "Ich war hier! !{workflow.projectDir}" > out.daa
   #diamond !{params.steps.annotation.diamond.mode} !{params.steps.annotation.diamond.cores} !{params.steps.annotation.diamond.sensitivity} !{params.steps.annotation.diamond.outfmt} --out diamond.!{sample}.daa --db databases/kegg-2021-01.dmnd --query input_diamond/*.fa !{params.steps.annotation.diamond.targetseq} !{params.steps.annotation.diamond.evalue
   '''
}

process pTest {
   // Difficult part. Without Docker, tail -1 should be used. With Docker the script changes and head -1 in addition to sed musst be used to remove an omnious \ bevore $PATH.
   beforeScript "eval \$(grep PATH .command.run | head -1 | sed 's/\\//g');eval \$(grep AWS_ACCESS .command.run | head -1);eval \$(grep AWS_SECRET .command.run | head -1);s5cmd ${params.steps.annotation.diamond.s5cmd} cp -u -s ${params.steps.annotation.diamond.database} ${params.databases}"
   
   // For some reason ${params.diamond_tag} does not get the value from the nextflow.config, despite it working with params.databases in the beforeScript.
   container "92bdeebcf94f"
 
   tag "$sample"

   label 'large'

   publishDir params.output, saveAs: { filename -> getOutput("${sample}", params.runid, "blast", filename) }

   input:
   tuple val(sample), val(project)
   output:
   tuple val("${sample}"), path("*.out"), emit: test_results

   shell:
   '''
   diamond help > help.out
   '''
}


/**
*
* This entry point uses a file to search all existing s3 projects with different similarity search tools listed therein.
* The .tsv file has to have: DATASET PATH entries. 
*
**/
workflow wAnnotateS3File {

   take:
      projectTableFile
   main:
      projectTableFile | _wAnnotation
}


workflow _wAnnotation {

   take:
      projectTableFile
   main:
      projectTableFile | splitCsv(sep: '\t', header: true) | view() | pTest | set { projectTable }
   emit:
      projectTable

   
}
