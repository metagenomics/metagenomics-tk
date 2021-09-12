nextflow.enable.dsl=2

MODULE="annotation"
VERSION="0.1.0"

def getOutput(SAMPLE, RUNID, TOOL, filename){
    return SAMPLE + '/' + RUNID + '/' + MODULE + '/' + VERSION + '/' + TOOL + '/' + filename
}


process pDiamond {
   // Difficult part. Without Docker, tail -1 should be used. With Docker the script changes and head -1 in addition to sed musst be used to remove an omnious \ bevore $PATH.
   //beforeScript "eval \$(grep PATH .command.run | head -1 | sed 's/\\//g');eval \$(grep AWS_ACCESS .command.run | head -1);eval \$(grep AWS_SECRET .command.run | head -1);s5cmd ${params.steps.annotation.diamond.s5cmd} cp -u -s ${params.steps.annotation.diamond.database} ${params.databases}"
   
   // For some reason ${params.diamond_tag} does not get the value from the nextflow.config, despite it working with params.databases in the beforeScript.
   container "quay.io/biocontainers/diamond:${params.diamond_tag}"
 
   tag "$sample"

   label 'large'

   publishDir params.output, saveAs: { filename -> getOutput("${sample}", params.runid, "diamond", filename) }

   // Input not nice, data should be collected beforehand
   input:
   tuple val(sample), file(BinFile)
   output:
   tuple val("${sample}"), path("diamond.${sample}.out"), emit: diamond_results

   shell:
   '''
   s5cmd !{params.steps.annotation.s5cmd} cp -u -s !{params.steps.annotation.diamond.database} !{params.databases}
   diamond !{params.steps.annotation.diamond.mode} !{params.steps.annotation.diamond.cores} !{params.steps.annotation.diamond.sensitivity} !{params.steps.annotation.diamond.outfmt} --out diamond.!{sample}.out --db !{params.databases}kegg-2021-01.dmnd --query !{BinFile} !{params.steps.annotation.diamond.targetseq} !{params.steps.annotation.diamond.evalue} 
   '''
}


process pKEGGFromDiamond {

   tag "$sample"

   label 'small'

   container "pbelmann/python-env:${params.python_env_tag}"

   publishDir params.output, saveAs: { filename -> getOutput("${sample}", params.runid, "keggFromDiamond", filename) }

   input:
   tuple val(sample), file(diamond_result)

   output:
   tuple val("${sample}"), path("${sample}_kegg.tsv"), emit: kegg_diamond_results

   shell:
   '''
   s5cmd !{params.steps.annotation.s5cmd} cp -u -s !{params.steps.annotation.kegg.database} !{params.databases}kegg
   diamond2kegg.py !{diamond_result} !{params.databases}kegg !{sample}_kegg.tsv 
   '''
}


/**
*
* This entry point uses a file to grab all .fa files in the referenced directoryis for annotation.
* The .tsv file has to have: DATASET PATH entries. 
*
**/
workflow wAnnotateLocalFile {

   take:
      projectTableFile
   main:
      projectTableFile | splitCsv(sep: '\t', header: true) | map{ it -> [it.DATASET, file(it.PATH)] } | collectFile() | map{ it -> [it.name, it]} | view () | _wAnnotation
}


workflow _wAnnotation {

   take:
      bins
   main:
      bins | pDiamond | pKEGGFromDiamond |set { projectTable }
   emit:
      projectTable

   
}
