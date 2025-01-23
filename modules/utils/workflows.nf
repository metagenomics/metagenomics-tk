import java.util.regex.*;

/*
* This method either returns file path or url
*/
def getPath(f){
  return params.input.startsWith("s3://")? "s3:/" + f: f
}

workflow _wGetCheckm {
  take:
    selectedSRAMagAttributes
  main:
    // get Checkm results
    Pattern checkmPattern = Pattern.compile('.*/magAttributes/' + params.modules.magAttributes.version.major + '..*/.*/.*_checkm_generated.tsv$')
    selectedSRAMagAttributes | filter({ sra, path -> checkmPattern.matcher(path.toString()).matches()}) \
     | set { checkmFiles }

    checkmFiles | splitCsv(header: ["SAMPLE", "BIN_ID", "Marker lineage", "# genomes", "# markers", \
          "# marker sets", "0", "1", "2", "3", "4", "5+", "COMPLETENESS", "CONTAMINATION", "HETEROGENEITY"], sep: '\t') \
     | map { sra, bins -> bins} \
     | set { checkm }

    // get Checkm2 results
    Pattern checkm2Pattern = Pattern.compile('.*/magAttributes/' + params.modules.magAttributes.version.major + '..*/.*/.*_checkm2_generated.tsv$')
    selectedSRAMagAttributes | filter({ sra, path -> checkm2Pattern.matcher(path.toString()).matches()}) \
     | set { checkm2Files }

    checkm2Files | splitCsv(header: true, sep: '\t') \
     | map { sra, bins -> bins} \
     | set { checkm2 }

    // We allow only to execute checkm or checkm2 but not both
    checkm \
       | mix(checkm2) \
       | set {checkm}

    checkmFiles \
       | mix(checkm2Files) \
       | set {checkmFiles}

  emit:
    checkm
    checkmFiles
}

workflow _wGetAssemblyFiles {
  take:
    assemblyFiles
  main:
    // get Illumina assembly files
    Pattern assemblyIlluminaPattern = Pattern.compile('.*/assembly/' + params.modules.assembly.version.major + '..*/.*/.*_contigs.fa.gz$')
    assemblyFiles | filter({ sra, path -> assemblyIlluminaPattern.matcher(path.toString()).matches()}) \
     | map{ sra,f -> [sra, getPath(f)] } | set { illuminaAssembly }

    // get ONT assembly files
    Pattern assemblyONTPattern = Pattern.compile('.*/assemblyONT/' + params.modules.assemblyONT.version.major + '..*/.*/.*_contigs.fa.gz$')
    assemblyFiles | filter({ sra, path -> assemblyONTPattern.matcher(path.toString()).matches()}) \
      | map{ sra,f -> [sra, getPath(f)] } | set { ontAssembly }
  emit:
    illuminaAssembly    
    ontAssembly
}

workflow _wGetIlluminaBinningFiles {
  take:
     binningFiles
  main:
      // Figure out if binrefinement via magscot was used.
      // Example Binning related file of the binningFiles channel:
      // [test2, /vol/spool/peter/meta-omics-toolkit/output/test2/1/binning/0.5.0/magscot]
      FILE_PATH_IDX = 1
      IS_MAGSCOT = "magscot"
      IS_NOT_MAGSCOT = "no_magscot"
      binningFiles | filter( path -> path[FILE_PATH_IDX].endsWith(IS_MAGSCOT)) | map { it -> IS_MAGSCOT } \
        | unique | ifEmpty(IS_NOT_MAGSCOT) | set { isMagscot }

      // If magscot was used then get all files that are placed in the magscot folder. If magscot was not used then get all binning related files
      MAGSCOT_FLAG = 2
      BINNING_FILE_PATH = 1
      SAMPLE_IDX = 0
      Pattern magscotPattern = Pattern.compile('.*/binning/' + params.modules.binning.version.major + '..*/magscot/.*$')
      binningFiles | combine(isMagscot) \
        | filter( binningFile -> binningFile[MAGSCOT_FLAG] == IS_MAGSCOT ? magscotPattern.matcher(binningFile[BINNING_FILE_PATH].toString()).matches() : true ) \
        | map { binFile -> [binFile[SAMPLE_IDX], binFile[BINNING_FILE_PATH]] } \
        | set { filteredBinningFiles  }
  emit:
      filteredBinningFiles
}
