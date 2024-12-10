include { pDumpLogs } from '../utils/processes'



def getOutput(SAMPLE, RUNID, TOOL, filename){
    return SAMPLE + '/' + RUNID + '/' + params.modules.export.name + '/' + 
           params.modules.export.version.major + "." +
           params.modules.export.version.minor + "." +
           params.modules.export.version.patch + 
           '/' + TOOL + '/' + filename
}

process pEMGBAnnotatedContigs {

    container "${params.emgbAnnotatedContigs_image}"

    tag "Sample: $sample"

    publishDir params.output, mode: "${params.publishDirMode}", saveAs: { filename -> getOutput("${sample}", params.runid, "emgb", filename) }

    when params.steps.containsKey("export") && params.steps.export.containsKey("emgb")

    containerOptions  " --entrypoint='' "

    label 'tiny'

    input:
    tuple val(sample), path(contigs), path(mapping), path("bins/*")

    output:

    tuple val("${sample}"), path("${sample}.contigs.json.gz"), emit: json
    tuple env(FILE_ID), val("${output}"), val(params.LOG_LEVELS.INFO), file(".command.sh"), \
	file(".command.out"), file(".command.err"), file(".command.log"), emit: logs

    shell:
    template 'emgbAnnotatedContigs.sh'
}


process pEMGBAnnotatedBins {

    container "${params.emgbAnnotatedBins_image}"

    tag "Sample: $sample"

    publishDir params.output, mode: "${params.publishDirMode}", saveAs: { filename -> getOutput("${sample}", params.runid, "emgb", filename) }

    when params.steps.containsKey("export") && params.steps.export.containsKey("emgb")

    containerOptions  " --entrypoint='' "

    label 'tiny'

    input:
    tuple val(sample), path(checkm), path(gtdbtk), path("bins/*")

    output:
    tuple val("${sample}"), path("${sample}.bins.json.gz"), emit: json
    tuple env(FILE_ID), val("${output}"), val(params.LOG_LEVELS.INFO), file(".command.sh"), \
	file(".command.out"), file(".command.err"), file(".command.log"), emit: logs

    shell:
    template 'emgbAnnotatedBins.sh'
}


process pEMGBAnnotatedGenes {

    container "${params.emgbAnnotatedGenes_image}"

    tag "Sample: $sample"

    publishDir params.output, mode: "${params.publishDirMode}", saveAs: { filename -> getOutput("${sample}", params.runid, "emgb", filename) }

    secret { "${S3_EMGB_TITLES_ACCESS}"!="" ? ["S3_EMGB_TITLES_ACCESS", "S3_EMGB_TITLES_SECRET"] : [] } 

    secret { "${S3_EMGB_KEGG_ACCESS}"!="" ? ["S3_EMGB_KEGG_ACCESS", "S3_EMGB_KEGG_SECRET"] : [] } 

    containerOptions Utils.getDockerMount(params.steps?.export?.emgb?.titles?.database, params) + Utils.getDockerMount(params.steps?.export?.emgb?.kegg?.database, params) + Utils.getDockerNetwork() + " --entrypoint='' "

    when params.steps.containsKey("export") && params.steps.export.containsKey("emgb")

    memory { Utils.getMemoryResources(params.resources.medium, "${sample}", task.attempt, params.resources) }

    cpus { Utils.getCPUsResources(params.resources.medium, "${sample}", task.attempt, params.resources) }

    beforeScript Utils.getCreateDatabaseDirCommand("${params.polished.databases}")

    input:
    tuple val(sample), path("gff/*"), path("ffn/*"), path("faa/*"), path("taxonomy/*"), path("blastResult/*"), path("blastKeggResult/*"), path("bins/*")

    output:
    tuple val("${sample}"), path("${sample}.genes.json.gz"), emit: json
    tuple env(FILE_ID), val("${output}"), val(params.LOG_LEVELS.INFO), file(".command.sh"), \
	file(".command.out"), file(".command.err"), file(".command.log"), emit: logs

    shell:
    TITLES_S5CMD_PARAMS=params?.steps?.export?.emgb?.titles?.database?.download?.s5cmd?.params ?: "" 
    TITLES_DOWNLOAD_LINK=params?.steps?.export?.emgb?.titles?.database?.download?.source ?: ""
    TITLES_MD5SUM=params.steps?.export?.emgb?.titles?.database?.download?.md5sum ?: ""
    TITLES_EXTRACTED_DB=params.steps?.export?.emgb?.titles?.database?.extractedDBPath ?: ""
    KEGG_S5CMD_PARAMS=params?.steps?.export?.emgb?.kegg?.database?.download?.s5cmd?.params ?: "" 
    KEGG_DOWNLOAD_LINK=params?.steps?.export?.emgb?.kegg?.database?.download?.source ?: ""
    KEGG_MD5SUM=params.steps?.export?.emgb?.kegg?.database?.download?.md5sum ?: ""
    KEGG_EXTRACTED_DB=params.steps?.export?.emgb?.kegg?.database?.extractedDBPath ?: ""
 
    TITLES_S3_EMGB_ACCESS=params?.steps?.export?.emgb?.titles?.database?.download?.s5cmd && TITLES_S5CMD_PARAMS.indexOf("--no-sign-request") == -1 ? "\$S3_EMGB_TITLES_ACCESS" : ""
    TITLES_S3_EMGB_SECRET=params?.steps?.export?.emgb?.titles?.database?.download?.s5cmd && TITLES_S5CMD_PARAMS.indexOf("--no-sign-request") == -1 ? "\$S3_EMGB_TITLES_SECRET" : ""

    KEGG_S3_EMGB_ACCESS=params?.steps?.export?.emgb?.kegg?.database?.download?.s5cmd && KEGG_S5CMD_PARAMS.indexOf("--no-sign-request") == -1 ? "\$S3_EMGB_KEGG_ACCESS" : ""
    KEGG_S3_EMGB_SECRET=params?.steps?.export?.emgb?.kegg?.database?.download?.s5cmd && KEGG_S5CMD_PARAMS.indexOf("--no-sign-request") == -1 ? "\$S3_EMGB_KEGG_SECRET" : ""
    template 'emgbAnnotatedGenes.sh'
}



// CONTIGS: [test1, /vol/spool/metagenomics-tk/work_wFullPipeline/cd/7c943c0808e6d36c72a64834e7a88e/test1_contigs.fa.gz]
// [test2, /vol/spool/metagenomics-tk/work_wFullPipeline/00/71fd0e590f4f03fda4cb87903a3b38/test2_contigs.fa.gz]
// mapping: [test2, /vol/spool/metagenomics-tk/work_wFullPipeline/7c/43a3513d49143a76a29c4404417997/test2.bam]
// bins: [test1, [/vol/spool/metagenomics-tk/work_wFullPipeline/95/085ffbe17a6f2ebfb60709d3f33cf3/test1_bin.1.fa, /vol/spool/metagenomics-tk/work_wFullPipeline/95/085ffbe17a6f2ebfb60709d3f33cf3/test1_bin.2.fa]]
// gtdbtk: [test2, [/vol/spool/metagenomics-tk/work_wFullPipeline/d2/5b5937c7a950d87cd796e2bf80cfcd/chunk_00_test2_gtdbtk.bac120.summary.tsv]]
// checkm: [test2, /vol/spool/metagenomics-tk/work_wFullPipeline/ff/1292ed3399aada6fa1dffc18cfbc12/test2_checkm2_EfY2cPIh.tsv]
// gff,ffn,faa: [test2, test2_bin.2.fa, /vol/spool/metagenomics-tk/work_wFullPipeline/78/1224783c7ae2d9e0a8d4264893e47a/test2_bin.2.gff.gz]
// mmseqsTaxonomy: [ncbi_nr, test1, /vol/spool/metagenomics-tk/work_wFullPipeline/93/56bf5c34d7524930d73a6ef45ae23d/test1_binned.ncbi_nr.taxonomy.tsv]
//mmseqsBlast:
// [bacmet20_predicted, test2, binned, 1, 2922, 1, /vol/spool/metagenomics-tk/work_wFullPipeline/3c/e5221be7fd305dc056356201de7d61/test2_binned.1.2922.bacmet20_predicted.blast.tsv]
// [uniref90, test1, binned, 1, 2923, 1, /vol/spool/metagenomics-tk/work_wFullPipeline/a0/8e7a2a7b928731dbcfd62cbf546d4f/test1_binned.1.2923.uniref90.blast.tsv]
// [bacmet20_predicted, test1, binned, 1, 2923, 1, /vol/spool/metagenomics-tk/work_wFullPipeline/2f/dffc66ee85b2a1baf296323fcc894c/test1_binned.1.2923.bacmet20_predicted.blast.tsv]
// [kegg, test2, binned, 1, 2922, 1, /vol/spool/metagenomics-tk/work_wFullPipeline/5a/74adca85c1b4005aa490440f8b05ba/test2_binned.1.2922.kegg.blast.tsv]
// [uniref90, test2, binned, 1, 2922, 1, /vol/spool/metagenomics-tk/work_wFullPipeline/57/5f4e7f5ae75937613d7c813be6c8fd/test2_binned.1.2922.uniref90.blast.tsv]
// [kegg, test1, binned, 1, 2923, 1, /vol/spool/metagenomics-tk/work_wFullPipeline/6b/41a5ce820f9fb5e897b6a9109e5678/test1_binned.1.2923.kegg.blast.tsv]
workflow wEMGBList {
  take:
    contigs
    mapping
    bins
    gtdbtk
    checkm
    gff
    ffn
    faa
    mmseqsTaxonomy
    mmseqsBlast
  main:
    SAMPLE_IDX = 0
    contigs |  combine(mapping, by: SAMPLE_IDX) | combine(bins, by: SAMPLE_IDX) | pEMGBAnnotatedContigs
    checkm | combine(gtdbtk, by: SAMPLE_IDX) | combine(bins, by: SAMPLE_IDX) | pEMGBAnnotatedBins

    mmseqsBlast | filter { db, sample, type, start, end, chunkNumber, blastResult -> db == "uniref90" } \
	| map { db, sample, type, start, end, chunkNumber, blastResult -> [sample, blastResult] } \
	| groupTuple(by: SAMPLE_IDX) | set { selectedDBBlastResults }

    mmseqsBlast | filter { db, sample, type, start, end, chunkNumber, blastResult -> db == "kegg" } \
	| map { db, sample, type, start, end, chunkNumber, blastResult -> [sample, blastResult] } \
	| groupTuple(by: SAMPLE_IDX) | set { selectedKeggBlastResults }

    gff | map { sample, bin, gff -> [sample, gff] } | groupTuple(by: SAMPLE_IDX) \
	| combine(ffn | map { sample, bin, ffn -> [sample, ffn] } | groupTuple(by: SAMPLE_IDX), by: SAMPLE_IDX) \
	| combine(faa | map { sample, bin, faa -> [sample, faa] } | groupTuple(by: SAMPLE_IDX), by: SAMPLE_IDX) \
	| combine(mmseqsTaxonomy | map { db, sample, blastResult -> [sample, blastResult] } | groupTuple(by: SAMPLE_IDX), by: SAMPLE_IDX) \
	| combine(selectedDBBlastResults, by: SAMPLE_IDX) \
	| combine(selectedKeggBlastResults, by: SAMPLE_IDX) \
	| combine(bins, by: SAMPLE_IDX) | pEMGBAnnotatedGenes
}
