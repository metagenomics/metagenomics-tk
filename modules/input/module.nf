nextflow.enable.dsl=2
import nextflow.splitter.CsvSplitter
import java.util.regex.*;


def getOutput(SAMPLE, RUNID, TOOL, filename){
    return SAMPLE + '/' + RUNID + '/' + params.modules.qc.name + '/' + 
           params.modules.qc.version.major + "." +
           params.modules.qc.version.minor + "." +
           params.modules.qc.version.patch + 
           '/' + TOOL + '/' + filename
}


process pGetSRAPath {

    errorStrategy 'retry'
    
    executor 'local'

    tag "Accession: $sraid"

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


process pGetMetadata {

    label 'tiny'

    container "${params.pysradb_image}"

    errorStrategy 'retry'

    time '3m'

    tag "Accession: $sraid"

    maxForks 3
    
    when:
    params?.input.containsKey("SRA")

    input:
    val(sraid)

    output:
    tuple val("${sraid}"), env(INSTRUMENT)

    shell:
    '''
    pysradb metadata !{sraid} --saveto output.tsv
    INSTRUMENT=$(tail -n 1  output.tsv | cut -f 17)
    '''
   
}

process pGetSRAIDs {

    label 'tiny'

    container "${params.pysradb_image}"

    errorStrategy 'retry'

    time '3m'

    maxForks 4

    tag "Accession: $sraid"
    
    when:
    params?.input.containsKey("SRA")

    input:
    val(sraid)

    output:
    path("output.tsv"), optional: true, emit: sraRunIDs
    env(NOT_FOUND_ID), optional: true, emit: notFoundID

    shell:
    '''
    pysradb metadata !{sraid} --saveto unfiltered_output.tsv
    TAB=`echo -e "\t"`
    cut -f 1,20 unfiltered_output.tsv | grep -e "$(printf '\t')!{sraid}$" -e "^!{sraid}$(printf '\t')" \
	| sed -r '/^\s*$/d'  > containsIDTest.txt
    if [ -s containsIDTest.txt ]; then
        echo "run_accession" > output.tsv
	cut -f 2 containsIDTest.txt >> output.tsv
    else 
        NOT_FOUND_ID=!{sraid}
	echo "No result for ID !{sraid} found.";
    fi
    '''
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
         BUCKET = params.input.SRA.S3.bucket

         Channel.fromPath(params.input.SRA.S3.path) \
		| set { idsFromPath }

         if(params.input.SRA.S3.watch){
             Channel.watchPath(params.input.SRA.S3.path, 'create,modify') \
		| set { idsFromWatch }
             idsFromWatch | mix(idsFromPath) | set {files}
         } else {
           idsFromPath | set { files }
         }

         files | splitCsv(sep: "\t", header: true) | map { it -> it.ACCESSION} | flatten | unique \
	       | pGetSRAIDs 

         pGetSRAIDs.out.sraRunIDs | splitCsv(sep: "\t", header: true) \
	        | map { sample -> sample.run_accession }  | unique \
		| branch {
                   passed: it.length() <= MAX_LENGTH && it.length() >= MIN_LENGTH
                           return it
                   failed: it.length() > MAX_LENGTH || it.length() < MIN_LENGTH
                           return it
	        } | set { filteredIDs }

         filteredIDs.passed | pGetSRAPath 

         pGetSRAPath.out.passed | map { it -> [ it[SAMPLE_IDX], file(BUCKET + it[FASTQ_FILES_IDX]).listFiles()]} \
                | map { it -> [ it[SAMPLE_IDX],  it[FASTQ_FILES_IDX].collect({ "s3:/$it" }) ] }
                | _wCheckSRAFiles

        emit:
          fastqs = _wCheckSRAFiles.out.passedSamples
          failedSRAFastqFiles = _wCheckSRAFiles.out.failedSRAIDs
          incorrectAccessions = filteredIDs.failed
          notFoundAccessions = pGetSRAIDs.out.notFoundID
}


/*
*
* Method to parse a tsv file containing a column with an ACCESSION header.
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
       runAccessions.add( row['ACCESSION'] )
    }
    return runAccessions
}


workflow _wCheckSRAFiles {
   take:
     samples
   main:
     SAMPLE_IDX = 0
     INSTRUMENT_IDX=1
     FASTQ_LEFT_IDX=2
     FASTQ_RIGHT_IDX=3
     FASTQ_FILES_IDX = 2
     ONT_FASTQ_FILE_IDX =2

     samples | map({sample -> sample[SAMPLE_IDX]}) | pGetMetadata \
        | combine(samples, by: SAMPLE_IDX) | branch {
         ONT: it[INSTRUMENT_IDX]=="OXFORD_NANOPORE"
         ILLUMINA: it[INSTRUMENT_IDX]=="ILLUMINA"
     } | set { samplesType }

     // Ensure that paired end files are included
       Pattern illuminaPattern = Pattern.compile(params.input?.SRA?.pattern.illumina);
       samplesType.ILLUMINA | branch {
            passed: it[FASTQ_FILES_IDX].size() >= 2 \
			&& it[FASTQ_FILES_IDX].stream().filter({ sraFile -> illuminaPattern.matcher(sraFile.normalize().toString()).matches() }).count() >= 2
                    return it
            other: true
                    return it
     } | set { sraFilesIllumina  }

     Pattern ontPattern = Pattern.compile(params.input?.SRA?.pattern.ont);
     samplesType.ONT | branch {
            passed: it[FASTQ_FILES_IDX].size() >= 1 \
			&& it[FASTQ_FILES_IDX].stream().anyMatch { sraFile -> ontPattern.matcher(sraFile.toString()).matches() }
                    return it
            other: true
                    return it
     } | set { sraFilesOnt  }

     // Ignore _3.fastq.gz files 
     sraFilesIllumina.passed | map({ sample ->  [sample[SAMPLE_IDX], sample[INSTRUMENT_IDX], \
	sample[FASTQ_FILES_IDX].findAll { sraFile -> illuminaPattern.matcher(sraFile.normalize().toString()).matches() } ].flatten()} ) \

     //Set map entries
     | map({ sample -> [SAMPLE:sample[SAMPLE_IDX], TYPE:sample[INSTRUMENT_IDX], READS1:sample[FASTQ_LEFT_IDX], READS2:sample[FASTQ_RIGHT_IDX]]}) \
     | set {illuminaFastqs}

      // Ignore _3.fastq.gz and _2.fastq.gz files 
      sraFilesOnt.passed | map({ sample ->  [sample[SAMPLE_IDX], sample[INSTRUMENT_IDX], \
	sample[FASTQ_FILES_IDX].findAll { sraFile -> ontPattern.matcher(sraFile.toString()).matches() } ].flatten()} ) \

     //Set map entries
     | map({ sample -> [SAMPLE:sample[SAMPLE_IDX], TYPE:sample[INSTRUMENT_IDX], READS:sample[ONT_FASTQ_FILE_IDX]]}) \
     | set {ontFastqs}

     illuminaFastqs | mix(ontFastqs) | set { fastqs }
     
     // Return ids of samples without two fastq files
     sraFilesIllumina.other | mix(sraFilesOnt.other) | map { it -> it[SAMPLE_IDX] } | set { failedSRAIDs }

    emit:
      passedSamples = fastqs
      failedSRAIDs = failedSRAIDs
}



workflow _wOntReads {
       main:
         Channel.fromPath(params.input.ont.path) \
                | set { idsFromPath }

         if(params.input.ont.watch){
            Channel.watchPath(params.input.ont.path, 'create,modify') \
                | set {idsFromWatch}

            idsFromWatch | mix(idsFromPath) | set { files }
         } else {
            idsFromPath | set {files}
         }

         files |  splitCsv(sep: '\t', header: true)  | unique | set { fastqs }
       emit:
         fastqs
}

workflow _wSRANCBI {
       main:
         MAX_LENGTH=12
         MIN_LENGTH=9

         // Parse TSV file to get access numbers
         accessions = fetchRunAccessions(params.input.SRA.NCBI.path)

         // check if the number of SRA files is correct and return the correct format 
         ACCESSION_ID = 0
         FASTQ_LIST = 1 
         Channel.fromSRA(accessions.unique()) | set { foundSRAFiles }
         foundSRAFiles | map { f -> [f[ACCESSION_ID],  Utils.asList(f[FASTQ_LIST])] } | _wCheckSRAFiles 
         Channel.from(accessions) \
		| combine(foundSRAFiles | map { sample -> sample[ACCESSION_ID]} | toList() | toList()) \
  		| filter({ id,idList -> !idList.contains(id) }) | map{ id -> id[ACCESSION_ID] } \
		| set { notFoundAccessions }

       emit:
         fastqs = _wCheckSRAFiles.out.passedSamples
         failedSRAFastqFiles = _wCheckSRAFiles.out.failedSRAIDs
         notFoundAccessions = notFoundAccessions
}


/*
 *  The input modules defined three input sources: SRA NCBI, generic SRA S3 source that contains a column consisting of SRA IDs 
 *  and a generic source that allows to consume a file containing local path, https or S3 links of paired end or nanopore reads.
 *  SRA NCBI and generic SRA S3 source must contain a column with `ACCESSION` column header and the file for the generic source 
 *  must contain the columns headers `SAMPLE`, `READS1` and `READS2` for paired end and `SAMPLE` and `READS` for nanopore data.
 * 
 *  In all cases a channel is returned containing values of the format: [TYPE: illumina or ont, SAMPLE:name of the sample, READS1: left read, READS2: right read],
 *  [TYPE: illumina or ont, SAMPLE:name of the sample, READS1: left read, READS2: right read] 
 */
workflow wInputFile {
  main:
    inputTypes = params.input.keySet()
    fastqs = Channel.empty()
    
    if("SRA" in inputTypes ){
        SRA_MODE=params.input.SRA.keySet()
        datasetsS3Fastqs = Channel.empty()
        datasetsNCBIFastqs = Channel.empty()
        datasetsS3IncorrectAccessions = Channel.empty()
        datasetsS3FailedSRAFastqFiles = Channel.empty()
        datasetsNCBIIncorrectAccessions = Channel.empty()
        datasetsNCBIFailedSRAFastqFiles = Channel.empty()
        datasetsS3notFoundAccessions = Channel.empty()
        datasetsNCBInotFoundAccessions = Channel.empty()

        if("S3" in SRA_MODE){
           _wSRAS3() | set { datasetsS3 }
           datasetsS3.fastqs | set { datasetsS3Fastqs }
           datasetsS3.incorrectAccessions | set { datasetsS3IncorrectAccessions }
           datasetsS3.failedSRAFastqFiles | set { datasetsS3FailedSRAFastqFiles }
           datasetsS3.notFoundAccessions | set { datasetsS3notFoundAccessions }
        } 

        if ("NCBI" in SRA_MODE){
           _wSRANCBI() | set { datasetsNCBI }
           datasetsNCBI.fastqs | set { datasetsNCBIFastqs }
           datasetsNCBI.failedSRAFastqFiles | set { datasetsNCBIFailedSRAFastqFiles }
           datasetsNCBI.notFoundAccessions | set { datasetsNCBInotFoundAccessions }
        }

        datasetsS3Fastqs | mix(datasetsNCBIFastqs) | set { sraFastqs }

        fastqs | mix(sraFastqs) | set { fastqs }

        datasetsS3IncorrectAccessions | set {incorrectAccessions}

        datasetsS3FailedSRAFastqFiles | mix(datasetsNCBIFailedSRAFastqFiles) | set {failedSRAFastqFiles}

	datasetsNCBInotFoundAccessions | mix(datasetsS3notFoundAccessions) | set {notFoundFiles}

        incorrectAccessions \
         | collectFile(newLine: true, seed: "ACCESSION", name: 'incorrectAccessions.tsv', storeDir: params.logDir)

        incorrectAccessions \
         | view { id -> "WARNING: The length of the ID $id is not between 9 and 12 characters and therefore will be skipped" }

        failedSRAFastqFiles \
         | collectFile(newLine: true, seed: "ACCESSIN", name: 'accessionsMissingFiles.tsv', storeDir: params.logDir)

        failedSRAFastqFiles \
         | view { id -> "WARNING: The SRA ID $id does neither provide illumina nor nanopore reads according to the specified pattern and therefore will be skipped" }

        notFoundFiles \
         | view { id -> "WARNING: The following ids could not be found: " + id }

    }

    if("paired" in inputTypes){
        _wSplitReads().fastqs \
                | map { data -> data["TYPE"] = "ILLUMINA"; data } \
                | set { pairedChannel }
        fastqs | mix(pairedChannel) | set { fastqs }
    }

    if("ont" in inputTypes) {
        _wOntReads().fastqs \
                | map { data -> data["TYPE"] = "OXFORD_NANOPORE"; data } \
                | set { ontChannel }

        fastqs | mix(ontChannel) | set { fastqs }
    }
  emit:
    data = fastqs
}
