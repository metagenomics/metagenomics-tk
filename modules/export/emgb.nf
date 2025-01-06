import java.util.regex.*;

include { pDumpLogs } from '../utils/processes'

include { wSaveSettingsList } from '../config/module'
include { _wGetCheckm; _wGetAssemblyFiles; _wGetIlluminaBinningFiles } from '../utils/workflows'
include { collectModuleFiles } from '../utils/processes'


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
    output = getOutput("${sample}", params.runid, "emgb", "")
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
    output = getOutput("${sample}", params.runid, "emgb", "")
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
 
    TITLES_S3_EMGB_ACCESS=params?.steps?.export?.emgb?.titles?.database?.download?.s5cmd \
	&& TITLES_S5CMD_PARAMS.indexOf("--no-sign-request") == -1 ? "\$S3_EMGB_TITLES_ACCESS" : ""
    TITLES_S3_EMGB_SECRET=params?.steps?.export?.emgb?.titles?.database?.download?.s5cmd \
	&& TITLES_S5CMD_PARAMS.indexOf("--no-sign-request") == -1 ? "\$S3_EMGB_TITLES_SECRET" : ""

    KEGG_S3_EMGB_ACCESS=params?.steps?.export?.emgb?.kegg?.database?.download?.s5cmd && KEGG_S5CMD_PARAMS.indexOf("--no-sign-request") == -1 ? "\$S3_EMGB_KEGG_ACCESS" : ""
    KEGG_S3_EMGB_SECRET=params?.steps?.export?.emgb?.kegg?.database?.download?.s5cmd && KEGG_S5CMD_PARAMS.indexOf("--no-sign-request") == -1 ? "\$S3_EMGB_KEGG_SECRET" : ""
    output = getOutput("${sample}", params.runid, "emgb", "")
    template 'emgbAnnotatedGenes.sh'
}

/*
* This workflow entry point allows to aggregate information of different samples.
* It will perform analysis steps such as dereplication, read mapping and co-occurrence.
* The input files are automatically fetched as long as they adhere to the pipeline specification document (see documentation).
*/
workflow _wExportPipeline {
    def input = params.input
    def runID = params.runid

    SAMPLE_IDX = 0

    // List all available SRAIDs
    Channel.from(file(input).list()) | filter({ path -> !(path ==~ /.*summary$/) && !(path ==~ /null$/) }) \
     | filter({ path -> !(path ==~ /.*AGGREGATED$/)}) \
     | set { sraDatasets }

    // Save config File
    wSaveSettingsList(sraDatasets)

    sraDatasets | map { sra ->  [sra, input + "/" + sra + "/" + runID + "/" ]} \
     | set {sraIDs}

    // List all files in sample directories
    sraIDs | flatMap { sraID, path -> collectModuleFiles(path, sraID, [params.modules.qc])} | set { qcFiles }
    sraIDs | flatMap { sraID, path -> collectModuleFiles(path, sraID, [params.modules.binning]) } \
	| set { binningFilesIllumina } 

    sraIDs | flatMap { sraID, path -> collectModuleFiles(path, sraID, [params.modules.binningONT])} \
        | set { binningONTFiles }

    binningFilesIllumina | mix(binningONTFiles) | set {binningFiles}

    sraIDs | flatMap { sraID, path -> collectModuleFiles(path, sraID, [params.modules.magAttributes])} \
	| set { selectedSRAMagAttributes}

    // Checkm files
    selectedSRAMagAttributes | _wGetCheckm 
    _wGetCheckm.out.checkmFiles | set { checkmFiles }

    sraIDs | flatMap { sraID, path -> collectModuleFiles(path, sraID, [params.modules.assemblyONT]) } | set { assemblyONTFiles } 
    sraIDs | flatMap { sraID, path -> collectModuleFiles(path, sraID, [params.modules.assembly]) } | set { assemblyIlluminaFiles }
    sraIDs | flatMap { sraID, path -> collectModuleFiles(path, sraID, [params.modules.annotation])} | set { selectedAnnotation}

    assemblyONTFiles | mix(assemblyIlluminaFiles) | _wGetAssemblyFiles 
    _wGetAssemblyFiles.out.illuminaAssembly | mix(_wGetAssemblyFiles.out.ontAssembly) | set { assembly } 

    // get Bins
    Pattern binsIlluminaPattern = Pattern.compile('.*/binning/' + params.modules.binning.version.major + '..*/.*/.*_bin.*.fa$')
    binningFiles | _wGetIlluminaBinningFiles | filter({ sra, path -> binsIlluminaPattern.matcher(path.toString()).matches()}) \
     | set{ illuminaBins }

    Pattern binsONTPattern = Pattern.compile('.*/binningONT/' + params.modules.binningONT.version.major + '..*/.*/.*_bin.*.fa$')
    binningFiles | filter({ sra, path -> binsONTPattern.matcher(path.toString()).matches()}) \
     | mix(illuminaBins) | groupTuple(by: SAMPLE_IDX) | set{ binFiles }

    MAX_BIN_COUNTER = 1000000
    // get not Binned gff files
    Pattern annotationNotBinnedGffPattern = Pattern.compile('.*/annotation/' + params.modules.annotation.version.major + '..*/prokka/.*_notBinned.gff.gz$' )
    selectedAnnotation | filter({ sra, path -> annotationNotBinnedGffPattern.matcher(path.toString()).matches()}) \
     | set { notBinnedGff }

    // get Bin gff files
    Pattern annotationGffPattern = Pattern.compile('.*/annotation/' + params.modules.annotation.version.major + '..*/prokka/.*_bin..*.gff.gz$' )
    selectedAnnotation | filter({ sra, path -> annotationGffPattern.matcher(path.toString()).matches()}) \
     | mix(notBinnedGff) | map { sra, path -> [sra, file(path).baseName, path, MAX_BIN_COUNTER] }  | set { gff }

    // get not binned ffn files
    Pattern annotationNotBinnedFfnPattern = Pattern.compile('.*/annotation/' + params.modules.annotation.version.major + '..*/prokka/.*_notBinned.ffn.gz$' )
    selectedAnnotation | filter({ sra, path -> annotationNotBinnedFfnPattern.matcher(path.toString()).matches()}) \
     | set { notBinnedFfn }

    // get ffn files
    Pattern annotationFfnPattern = Pattern.compile('.*/annotation/' + params.modules.annotation.version.major + '..*/prokka/.*_bin..*.ffn.gz$' )
    selectedAnnotation | filter({ sra, path -> annotationFfnPattern.matcher(path.toString()).matches()}) \
     | mix(notBinnedFfn) | map { sra, path -> [sra, file(path).baseName, path, MAX_BIN_COUNTER] } | set { ffn }

    // not Binned faa files
    Pattern annotationNotBinnedPattern = Pattern.compile('.*/annotation/' + params.modules.annotation.version.major + '..*/prokka/.*_notBinned.faa.gz$' )
    selectedAnnotation | filter({ sra, path -> annotationNotBinnedPattern.matcher(path.toString()).matches()}) \
     | set { notBinnedFaa }

    // get Bin faa files
    Pattern annotationFnaPattern = Pattern.compile('.*/annotation/' + params.modules.annotation.version.major + '..*/prokka/.*_bin..*.faa.gz$' )
    selectedAnnotation | filter({ sra, path -> annotationFnaPattern.matcher(path.toString()).matches()}) \
     | mix(notBinnedFaa) | map { sra, path -> [sra, file(path).baseName, path, MAX_BIN_COUNTER] } | set { faa }

    def TAXONOMY_DB = "gtdb"
    def BLAST_DB = "uniref90"
    if(params.steps.containsKey("export") \
	&& params.steps.export.containsKey("emgb") \
        && params.steps.export.emgb.containsKey("additionalParams")){

	if(!params.steps.export.emgb.additionalParams.taxonomyDB.isEmpty()){
		TAXONOMY_DB = params.steps.export.emgb.additionalParams.taxonomyDB
	}

	if(!params.steps.export.emgb.additionalParams.blastDB.isEmpty()){
		BLAST_DB = params.steps.export.emgb.additionalParams.blastDB
	}
    }

    // get MMseqs unbinned taxonomy files
    Pattern annotationMmseqsUnbinnedTaxonomyPattern = Pattern.compile('.*/annotation/' + params.modules.annotation.version.major + \
	'..*/mmseqs2_taxonomy/.*/.*_unbinned.' + TAXONOMY_DB + '.taxonomy.tsv$' )
    selectedAnnotation | filter({ sra, path -> annotationMmseqsUnbinnedTaxonomyPattern.matcher(path.toString()).matches()}) \
     | map { sra, path -> [sra, TAXONOMY_DB, path, MAX_BIN_COUNTER]  } | set { mmseqsUnbinnedTaxonomy }

    // get MMseqs binned taxonomy files   output/test2/1/annotation/1.0.0/mmseqs2/ncbi_nr/   test2_unbinned.ncbi_nr.105001.112000.blast.tsv
    Pattern annotationMmseqsBinnedTaxonomyPattern = Pattern.compile('.*/annotation/' + params.modules.annotation.version.major \
	+ '..*/mmseqs2_taxonomy/.*/.*_binned.' + TAXONOMY_DB + '.taxonomy.tsv$' )
    selectedAnnotation | filter({ sra, path -> annotationMmseqsBinnedTaxonomyPattern.matcher(path.toString()).matches()}) \
     | map { sra, path -> [sra, TAXONOMY_DB, path, MAX_BIN_COUNTER] } | mix(mmseqsUnbinnedTaxonomy) | set { mmseqsTaxonomy }

    // get MMseqs unbinned blast files
    Pattern annotationMmseqsUnirefUnbinnedPattern = Pattern.compile('.*/annotation/' + params.modules.annotation.version.major + '..*/mmseqs2/.*/.*_unbinned.' + BLAST_DB + '.blast.tsv$' )
    selectedAnnotation | filter({ sra, path -> annotationMmseqsUnirefUnbinnedPattern.matcher(path.toString()).matches()}) \
     | map { sra, path -> [sra, BLAST_DB, path, MAX_BIN_COUNTER] } | set { mmseqsUnirefUnbinnedBlast }

    // get MMseqs binned blast files   output/test2/1/annotation/1.0.0/mmseqs2/ncbi_nr/   test2_unbinned.ncbi_nr.105001.112000.blast.tsv
    Pattern annotationMmseqsUnirefBinnedPattern = Pattern.compile('.*/annotation/' + params.modules.annotation.version.major + '..*/mmseqs2/.*/.*_binned.' + BLAST_DB + '.blast.tsv$' )
    selectedAnnotation | filter({ sra, path -> annotationMmseqsUnirefBinnedPattern.matcher(path.toString()).matches()}) \
     | map { sra, path -> [sra, BLAST_DB, path, MAX_BIN_COUNTER] } \
     | mix(mmseqsUnirefUnbinnedBlast) | set { mmseqsUnirefBlast }

    // get MMseqs unbinned kegg blast files
    Pattern annotationMmseqsKeggUnbinnedPattern = Pattern.compile('.*/annotation/' + params.modules.annotation.version.major + '..*/mmseqs2/.*/.*_unbinned.kegg.blast.tsv$' )
    selectedAnnotation | filter({ sra, path -> annotationMmseqsKeggUnbinnedPattern.matcher(path.toString()).matches()}) \
     | map { sra, path -> [sra, "kegg", path, MAX_BIN_COUNTER] } \
     | set { mmseqsKeggUnbinnedBlast }

    // get MMseqs binned blast files
    Pattern annotationMmseqsKeggBinnedPattern = Pattern.compile('.*/annotation/' + params.modules.annotation.version.major + '..*/mmseqs2/.*/.*_binned.kegg.blast.tsv$' )
    selectedAnnotation | filter({ sra, path -> annotationMmseqsKeggBinnedPattern.matcher(path.toString()).matches()}) \
     | map { sra, path -> [sra, "kegg", path, MAX_BIN_COUNTER] } \
     | mix(mmseqsKeggUnbinnedBlast) \
     | set { mmseqsKeggBlast }

    // get gtdbtk summary files
    Pattern gtdbSummaryPattern = Pattern.compile('.*/magAttributes/' + params.modules.magAttributes.version.major + '..*/.*/.*_gtdbtk_generated_summary_raw_combined.tsv$' )
    selectedSRAMagAttributes | filter({ sra, path -> gtdbSummaryPattern.matcher(path.toString()).matches()}) \
     | groupTuple(by: SAMPLE_IDX) | set { gtdbSummaryFiles }

    // get Mapping
    Pattern mappingIlluminaPattern = Pattern.compile('.*/binning/' + params.modules.binning.version.major + '..*/.*/.*.bam$')
    binningFiles | filter({ sra, path -> mappingIlluminaPattern.matcher(path.toString()).matches()}) \
     | set{ illuminaMapping }

    // get ONT mapping files
    Pattern mappingONTPattern = Pattern.compile('.*/binningONT/' + params.modules.binningONT.version.major + '..*/.*/.*.bam$')
    binningFiles | filter({ sra, path -> mappingONTPattern.matcher(path.toString()).matches()}) \
     | set{ ontMapping }

    illuminaMapping | mix(ontMapping) | set { mapping } 

    wEMGBList(assembly, \
	mapping, \
        binFiles, \
        gtdbSummaryFiles, \
        checkmFiles, \
        gff, \
        ffn, \
        faa, \
        mmseqsTaxonomy, \
        mmseqsUnirefBlast | mix(mmseqsKeggBlast) \
    )
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

    selectedDBBlastResults = Channel.empty()
    if(params.steps.containsKey("export") \
	&& params.steps.export.containsKey("emgb") \
        && params.steps.export.emgb.containsKey("additionalParams")){
          
        BLAST_DB = "ncbi_nr"
	if(!params.steps.export.emgb.additionalParams.blastDB.isEmpty()){
		BLAST_DB = params.steps.export.emgb.additionalParams.blastDB
	}
    	mmseqsBlast | filter { sample, db, blastResult, counter -> db == BLAST_DB } \
		| map { sample, db, blastResult, counter -> tuple(groupKey(sample.toString(), counter), blastResult) } \
		| groupTuple(remainder: true) | map { sample, blastResult -> [sample.toString(), blastResult] } | set { selectedDBBlastResults }
    }

    selectedTaxonomyDBResult = Channel.empty()
    if(params.steps.containsKey("export") \
	&& params.steps.export.containsKey("emgb") \
        && params.steps.export.emgb.containsKey("additionalParams")){
        TAXONOMY_DB = "gtdb"
	if(!params.steps.export.emgb.additionalParams.taxonomyDB.isEmpty()){
		TAXONOMY_DB = params.steps.export.emgb.additionalParams.taxonomyDB
	}
        mmseqsTaxonomy | filter { sample, db, blastRestults, counter -> db==TAXONOMY_DB } \
		| map { sample, db, blastResults, counter -> tuple(groupKey(sample.toString(), counter), blastResults) } | groupTuple(remainder: true) \
		| map { sample, blastResults -> [sample.toString(), blastResults]} | set { selectedTaxonomyDBResult }
    }

    mmseqsBlast | filter { sample, db, blastResult, counter -> db == "kegg" } \
		| map { sample, db, blastResult, counter -> tuple(groupKey(sample, counter), blastResult) } \
		| groupTuple(remainder: true) | map { sample, blastResult -> [sample.toString(), blastResult] } | set { selectedKeggBlastResults }

    gff | map { sample, bin, gff, counter -> tuple(groupKey(sample, counter), gff) } | groupTuple(remainder: true) \
	| combine(ffn | map { sample, bin, ffn, counter -> tuple(groupKey(sample, counter), ffn) } | groupTuple(remainder: true), by: SAMPLE_IDX) \
	| combine(faa | map { sample, bin, faa, counter -> tuple(groupKey(sample, counter), faa) } | groupTuple(remainder: true), by: SAMPLE_IDX) \
		| map { sample, ffn, faa, gff -> [sample.toString(), ffn, faa, gff] } \
	| combine(selectedTaxonomyDBResult, by: SAMPLE_IDX) \
	| combine(selectedDBBlastResults, by: SAMPLE_IDX) \
	| combine(selectedKeggBlastResults, by: SAMPLE_IDX) \
	| combine(bins, by: SAMPLE_IDX) | pEMGBAnnotatedGenes

    pEMGBAnnotatedContigs.out.logs | mix(pEMGBAnnotatedBins.out.logs, pEMGBAnnotatedGenes.out.logs) | pDumpLogs
}
