nextflow.enable.dsl=2
import nextflow.splitter.CsvSplitter


def getOutput(SAMPLE, RUNID, TOOL, filename){
    return SAMPLE + '/' + RUNID + '/' + params.modules.qc.name + '/' + 
           params.modules.qc.version.major + "." +
           params.modules.qc.version.minor + "." +
           params.modules.qc.version.patch + 
           '/' + TOOL + '/' + filename
}


process getSRAPath {

    errorStrategy 'retry'
    
    executor 'local'
    
    when:
    params?.input.containsKey("SRA") && params?.input.SRA.containsKey("S3") 

    input:
    val(sraid)

    output:
    tuple val("${sraid}"), val(path), optional: true, emit: passed

    exec:
    PREFIX=params.input.SRA.S3.prefix
    switch(sraid.length()) { 
      case 12:
        path=PREFIX + sraid.substring(0,6) + sraid.substring(9,12) + sraid
        break;
      case 11:
        path=PREFIX + sraid.substring(0,6) + "/0" + sraid.substring(9,11) + "/" + sraid
        break;
      case 10:
        path=PREFIX + sraid.substring(0,6) + "/00" + sraid.substring(9,10) + "/" + sraid
        break;
      case 9:
        path=PREFIX + sraid.substring(0,6) + "/" + sraid
        break;
    }   
}


workflow _wSplitReads {
       main:
         Channel.fromPath(params.input.paired.path) \
		| set { idsFromPath }

         if(params.input.paired.watch){
            Channel.watchPath(params.input.paired.path, 'create,modify') \
                | set {idsFromWatch}
            
            idsFromWatch | mix(idsFromPath) | set { files }
         } else {
            idsFromPath | set {files}
         } 

         files |  splitCsv(sep: '\t', header: true)  | unique | set { fastqs } 
       emit:
         fastqs
}



workflow _wSRAS3 {
       main:
         MAX_LENGTH=12
         MIN_LENGTH=9
         FASTQ_FILES_IDX = 1
         SAMPLE_IDX = 0
         BUCKET = "s3://ftp.era.ebi.ac.uk" 

         Channel.fromPath(params.input.SRA.S3.path) \
		| set { idsFromPath }

         if(params.input.SRA.S3.watch){
             Channel.watchPath(params.input.SRA.S3.path, 'create,modify') \
		| set { idsFromWatch }
             idsFromWatch | mix(idsFromPath) | set {files}
         } else {
           idsFromPath | set { files }
         }

         files | splitCsv(sep: "\t", header: true) | map { it -> it.RUN_ID} | flatten | unique  \
		| branch {
                   passed: it.length() <= MAX_LENGTH && it.length() >= MIN_LENGTH
                           return it
                   failed: it.length() > MAX_LENGTH || it.length() < MIN_LENGTH
                           return it
	        } | set { filteredIDs }

         filteredIDs.passed | getSRAPath 
 
         getSRAPath.out.passed | map { it -> [ it[SAMPLE_IDX], file(BUCKET + it[FASTQ_FILES_IDX]).listFiles()]} \
                | map { it -> [ it[SAMPLE_IDX],  it[FASTQ_FILES_IDX].collect({ "s3:/$it" }) ] }
                | _wCheckSRAFiles

        emit:
          fastqs = _wCheckSRAFiles.out.passedSamples
          failedSRAFastqFiles = _wCheckSRAFiles.out.failedSRAIDs
          incorrectAccessions = filteredIDs.failed
}


/*
*
* Method to parse a tsv file containing a column with an RUN_ID header.
* Input: Path to a TSV file
* Output: List of SRA run ids
*
*/
def fetchRunAccessions( tsv ) {

    def splitter = new CsvSplitter().options( header:true, sep:'\t' )
    def reader = new BufferedReader( new FileReader( tsv ) )

    splitter.parseHeader( reader )

    List<String> runAccessions = []
    Map<String, String> row

    while( row = splitter.fetchRecord( reader ) ) {
       runAccessions.add( row['RUN_ID'] )
    }
    return runAccessions
}


workflow _wCheckSRAFiles {
   take:
     samples
   main:
     SAMPLE_IDX = 0
     FASTQ_LEFT_IDX=1
     FASTQ_RIGHT_IDX=2
     FASTQ_FILES_IDX = 1

     // Ensure that at least two fastq files are in the folder
     samples | filter({sample -> sample[FASTQ_FILES_IDX].size() >= 2 }) \
         // Ensure that paired end files are included
	 | branch {
            passed: it[FASTQ_FILES_IDX].stream().allMatch { sraFile -> sraFile ==~ /.+(_1|_2).+$/ }
                    return it
            other: true
                    return it
         } | set { sraFiles  }
         // Ignore _3.fastq.gz files 
         sraFiles.passed | map({ sample ->  [sample[SAMPLE_IDX], sample[FASTQ_FILES_IDX].findAll { it ==~ /.+(_1|_2).+$/} ].flatten()} ) \
	 //Set map entries
	 | map({ sample -> [SAMPLE:sample[SAMPLE_IDX],READS1:sample[FASTQ_LEFT_IDX], READS2:sample[FASTQ_RIGHT_IDX]]}) \
	 | set {fastqs}

         // Return ids of samples without two fastq files
         sraFiles.other | map { it -> it[SAMPLE_IDX] } | set { failedSRAIDs }

         // report failed SRA ids
         failedSRAIDs \
	 | view { id -> "The following sample does not have two fastq files that match the patter '/.+(_1|_2).+\$/': $id " }

    emit:
      passedSamples = fastqs
      failedSRAIDs = failedSRAIDs
}



workflow _wSRANCBI {
       main:
         MAX_LENGTH=12
         MIN_LENGTH=9
          
         // Parse TSV file to get access numbers
         accessions = fetchRunAccessions(params.input.SRA.NCBI.path)

         // check that SRA IDs have the correct length
         accessionsFiltered = accessions.findAll{ id -> id.length() <= MAX_LENGTH && id.length() >= MIN_LENGTH }

         // get SRA IDs with incorrect length
         incorrectAccessions = Channel.from(accessions.findAll{ id -> id.length() > MAX_LENGTH || id.length() < MIN_LENGTH })

         // check if the number of SRA files is correct and return the correct format 
         Channel.fromSRA(accessionsFiltered.unique()) | _wCheckSRAFiles

       emit:
         fastqs = _wCheckSRAFiles.out.passedSamples
         failedSRAFastqFiles = _wCheckSRAFiles.out.failedSRAIDs
         incorrectAccessions = incorrectAccessions
}


/*
 *  The input modules defined three input sources: SRA NCBI, generic SRA S3 source that contains a column consisting of SRA IDs 
 *  and a generic source that allows to consume a file containing local path, https or S3 links of paired reads.
 *  SRA NCBI and generic SRA S3 source must contain a column with `RUN_ID` column header and the file for the generic source 
 *  must contain the columns headers `SAMPLE`, `READS1` and `READS2`.
 * 
 *  In all cases a channel is returned containing values of the format: [SAMPLE:name of the sample, READS1: left read, READS2: right read] 
 */
workflow wInputFile {
  main:
    switch(params.input.keySet()[0]){
      case "SRA":
        SRA_MODE=params.input.SRA.keySet()[0]
        if(SRA_MODE=="S3"){
           _wSRAS3() | set { datasets }
        } else if (SRA_MODE=="NCBI"){
           _wSRANCBI() | set { datasets }
        }
        break;
      case "paired":
        _wSplitReads() | set { datasets }
        break;
      default:
        println "WARNING: unknown or no input parameter specified!"
    }

    datasets.incorrectAccessions \
	| collectFile(newLine: true, seed: "RUN_ID", name: 'incorrectAccessions.tsv', storeDir: params.logDir)

    datasets.incorrectAccessions \
	| view { id -> "WARNING: The length of the ID $id is not between 9 and 12 characters and therefore will be skipped" }

    datasets.failedSRAFastqFiles \
	| collectFile(newLine: true, seed: "RUN_ID", name: 'accessionsMissingFiles.tsv', storeDir: params.logDir)

    datasets.failedSRAFastqFiles \
	| view { id -> "WARNING: The SRA ID $id does not provide paired end reads and therefore will be skipped" }

  emit:
    data = datasets.fastqs
}
