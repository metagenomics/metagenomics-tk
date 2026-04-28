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


process pGetSRAMetadataFromRemote {

    label 'tiny'

    container "${params.pysradb_image}"

    errorStrategy 'retry'

    time '3m'

    tag "Accession: $sraid"

    maxForks 3

    maxErrors 5
    
    when:
    params?.input.containsKey("SRA")

    input:
    val(sraid)

    output:
    tuple val("${sraid}"), env(INSTRUMENT)

    shell:
    '''
    pysradb metadata !{sraid} --saveto pysradb_output.tsv
    INSTRUMENT=$(tail -n 1 pysradb_output.tsv | cut -f 17)
    '''
}

/*
* If possible fetch from local SRA DB instead of fetching from NCBI.
* Input: List of SRA IDs
* Output: 
*  * foundSraRunIDs: File containing matadata of found SRA IDs.
*  * notFoundSraRunIDs: File containing SRA IDs of not found SRA IDs in DB.
*/
process pGetSRAIDsFromDB {

    label 'tiny'

    container "${params.ubuntu_image}"

    errorStrategy 'retry'

    when:
    params?.input.containsKey("SRA")

    input:
    val("sraids")

    output:
    path("found.csv"), emit: foundSraRunIDs
    path("not_found.csv"), emit: notFoundSraRunIDs

    shell:
    skipDB=(params.input.SRA.containsKey("skipDB") && params.input.SRA.skipDB) ? "true" : ""
    template("getSRAIDsFromDB.sh")
}

process pGetSRAIDsFromRemote {

    label 'tiny'

    container "${params.pysradb_image}"

    errorStrategy 'retry'

    time '10m'

    maxForks 4

    maxErrors 5

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
    TAB=`echo -e "\t"`
    count=1
    # Retry pysradb 3 times in case it can not find the correct id
    while [ $count -le 4 ]
    do
      pysradb metadata !{sraid} --saveto unfiltered_output.tsv
      cut -f 1,20 unfiltered_output.tsv | grep -e "$(printf '\t')!{sraid}$" -e "^!{sraid}$(printf '\t')" \
  	  | sed -r '/^\s*$/d'  > containsIDTest.txt

      if [ -s containsIDTest.txt ]; then
        echo "run_accession" > output.tsv
	    cut -f 2 containsIDTest.txt >> output.tsv
        unset NOT_FOUND_ID
        break;
      else 
	    echo "No result for ID !{sraid} found. Number of attempt: ${count}";
        sleep 2
        NOT_FOUND_ID=!{sraid}
        ((count++)) 
      fi
    done
    '''
}


/**
*
* This workflow processes Illumina input files provided via a sample sheet.
*
*/
workflow _wSplitReadsSheet {
       main:
         Channel.fromPath(params.input.paired.sheet) \
		    |  set { idsFromPath }

         if(params.input.paired.watch){
            Channel.watchPath(params.input.paired.sheet, 'create,modify') \
                | set {idsFromWatch}
            
            idsFromWatch | mix(idsFromPath) | set { files }
         } else {
            idsFromPath | set {files}
         } 
         files |  splitCsv(sep: '\t', header: true) | unique | set { fastqs } 
       emit:
         fastqs
}


/**
*
* This workflow processes Illumina input files provided via CLI.
*
*/
workflow _wSplitReadsFiles {
       main:
         def r1 = params.input.paired.r1.tokenize(' ')
         def r2 = params.input.paired.r2.tokenize(' ')
         def names = params.input.paired.names.tokenize(" ")

         if (r1.size() != r2.size() && r2.size() == r3.size() ) {
            error "Mismatch detected: --input.paired.r1, --input.paired.r2 and --input.paired.names should have the same number of values."
         }

         SAMPLE_IDX=0
         READS1_IDX=1
         READS2_IDX=2
         CO_BINNING_IDX=3

      

	     def samples = []
         fastqs = channel.empty()
         if(params.input.paired.containsKey("binGroup")){
            def binGroup = params.input.paired.binGroup.tokenize(" ")
            if (r1.size() != binGroup.size()) {
                error "Mismatch detected: --input.paired.r and --input.paired.binGroup should have the same number of values."
            }
	        samples = [names, r1, r2, binGroup].transpose()
            channel.from(samples) 
	            | map { sample -> [SAMPLE:sample[SAMPLE_IDX],READS1:file(sample[READS1_IDX]),READS2:file(sample[READS2_IDX]),CO_BINNING:sample[CO_BINNING_IDX]] }
                | set { fastqs }
         } else {
	        samples = [names, r1, r2].transpose()
            channel.from(samples) 
	            | map { sample -> [SAMPLE:sample[SAMPLE_IDX],READS1:file(sample[READS1_IDX]),READS2:file(sample[READS2_IDX])] }
                | set { fastqs }
         }


      emit:
         fastqs
}

workflow _wSRAS3 {
       main:
         MAX_LENGTH=12
         MIN_LENGTH=9
         FASTQ_FILES_IDX = 1
         SAMPLE_IDX = 0
         SAMPLE_CONTENT_IDX = 1
         CO_BINNING_IDX = 2
         BUCKET = params.input.SRA.S3.bucket

        files = Channel.empty()
        //If files are provided via a sample sheet we allow to watch the file.
        if("sheet" in params.input.SRA.S3){
            Channel.fromPath(params.input.SRA.S3.sheet) \
		    | set { idsFromPath }

            if(params.input.SRA.S3.watch){
                Channel.watchPath(params.input.SRA.S3.sheet, 'create,modify') \
		         | set { idsFromWatch }
                idsFromWatch | mix(idsFromPath) | set {files}
            } else {
                idsFromPath | set { files }
            }
        }

        // Files provided via CLI 
        idsFromCLIChannel = Channel.empty()
        def idsCLI = null
        if("id" in params.input.SRA.S3){
            idsCLI = params.input.SRA.S3.id.tokenize(" ")
            idsFromCLIChannel = Channel.from(params.input.SRA.S3.id.tokenize(" "))
        }

        binGroupFromCLIChannel = Channel.empty()
        if("binGroup" in params.input.SRA.S3){
            groupsCLI = params.input.SRA.S3.binGroup.tokenize(" ")
            if (idsCLI.size() != groupsCLI.size()) {
                error "Mismatch detected: --input.SRA.S3.id and --input.SRA.S3.binGroup should have the same number of values."
            }
	        coBinningSamplesCLI = [idsCLI, groupsCLI].transpose()
            binGroupFromCLIChannel = channel.from(coBinningSamplesCLI)
        }

        // First try to fetch SRA IDs from NCBI SRA DB
        BUFFER_SIZE_WATCH = 1 
        BUFFER_SIZE_DEFAULT = 50  
        files | splitCsv(sep: "\t", header: true) | map { it -> it.ACCESSION} 
            | mix(idsFromCLIChannel)
            | flatten | unique 
            | buffer(size: params.input.SRA.S3.watch ? BUFFER_SIZE_WATCH : BUFFER_SIZE_DEFAULT, remainder: true) 
            | pGetSRAIDsFromDB

        pGetSRAIDsFromDB.out.foundSraRunIDs 
            | splitCsv(sep: ",", header: true) 
            | map { sample -> sample.run_accession}
            | set { sraRunIDsFoundInDB }

        // SRA IDs not found in DB should be checked against NCBI directly
	    pGetSRAIDsFromDB.out.notFoundSraRunIDs 
	        | splitCsv(sep: ",", header: true) 
            | map { sample -> sample.run_accession}
	        | pGetSRAIDsFromRemote 

        pGetSRAIDsFromRemote.out.sraRunIDs | splitCsv(sep: "\t", header: true) \
 	        | map { sample -> sample.run_accession }  | unique \
            | set { sraIDsFoundRemote }


        // Mix both outputs
        sraIDsFoundRemote | mix(sraRunIDsFoundInDB) 
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

        files | splitCsv(sep: "\t", header: true) 
            | branch { sample ->
                coBinningSample: sample.containsKey("CO_BINNING") 
                notCoBinnedSample: !sample.containsKey("CO_BINNING") 
            } | set{ possibleCoBinnedSamples } 

         possibleCoBinnedSamples.coBinningSample
            | map { sample -> [sample.ACCESSION, sample.CO_BINNING] } 
            | mix(binGroupFromCLIChannel) | set { coBinningSamples }

        _wCheckSRAFiles.out.passedSamples | map { sample -> [sample.SAMPLE, sample]} 
            | combine(coBinningSamples, by: SAMPLE_IDX) 
            | map { sample -> sample[SAMPLE_CONTENT_IDX] + ["CO_BINNING":sample[CO_BINNING_IDX]] }   
            | mix(possibleCoBinnedSamples.notCoBinnedSample)
            | set { passedSamples }

        emit:
          fastqs = passedSamples
          failedSRAFastqFiles = _wCheckSRAFiles.out.failedSRAIDs
          incorrectAccessions = filteredIDs.failed
          notFoundAccessions = pGetSRAIDsFromRemote.out.notFoundID
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

    List runAccessions = []
    Map<String, String> row

    while( row = splitter.fetchRecord( reader ) ) {
       if(row.containsKey("CO_BINNING")){
        runAccessions.add( [row['ACCESSION'], row['CO_BINNING']] )
       } else {
        runAccessions.add([row['ACCESSION']])
       }
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

     // First try to fetch SRA IDs from NCBI SRA DB
     BUFFER_SIZE_WATCH = 1 
     BUFFER_SIZE_DEFAULT = 50  
     samples | map({sample -> sample[SAMPLE_IDX]}) 
        | buffer(size: params.input.SRA.S3.watch ? BUFFER_SIZE_WATCH : BUFFER_SIZE_DEFAULT, remainder: true) 
        | pGetSRAIDsFromDB

     pGetSRAIDsFromDB.out.foundSraRunIDs 
            | splitCsv(sep: ",", header: true) 
            | map { sample -> [ sample.run_accession, sample.instrument ]}
            | set { sraRunIDsMetadataFoundInDB }

     // SRA IDs not found in DB should be checked against NCBI directly
     pGetSRAIDsFromDB.out.notFoundSraRunIDs 
            | splitCsv(sep: ",", header: true) 
            | map { sample -> sample.run_accession}
            | pGetSRAMetadataFromRemote | set { sraRunIDsMetadataFoundRemote }


     // Mix both outputs
     sraRunIDsMetadataFoundInDB | mix(sraRunIDsMetadataFoundRemote)
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



workflow _wOntReadsSheet {
       main:
         Channel.fromPath(params.input.ont.sheet) \
                | set { idsFromSheet }

         if(params.input.ont.watch){
            Channel.watchPath(params.input.ont.sheet, 'create,modify') \
                | set {idsFromWatch}

            idsFromWatch | mix(idsFromSheet) | set { files }
         } else {
            idsFromSheet | set {files}
         }

         files |  splitCsv(sep: '\t', header: true)  | unique | set { fastqs }
       emit:
         fastqs
}

workflow _wOntReadsFiles {
       main:
         def r = params.input.ont.r.tokenize(' ')
         def names = params.input.ont.names.tokenize(" ")

         if (r.size() != names.size()) {
            error "Mismatch detected: --input.ont.r and --input.ont.names should have the same number of values."
         }

         SAMPLE_IDX=0
         READS_IDX=1
         CO_BINNING_IDX=2

	     def samples = []
         fastqs = channel.empty()
         if(params.input.ont.containsKey("binGroup")){
            def binGroup = params.input.ont.binGroup.tokenize(" ")
            if (r.size() != binGroup.size()) {
                error "Mismatch detected: --input.ont.r and --input.ont.binGroup should have the same number of values."
            }

	        samples = [names, r, binGroup].transpose()
            channel.from(samples) 
	            | map { sample -> [SAMPLE:sample[SAMPLE_IDX],READS:file(sample[READS_IDX]),CO_BINNING:sample[CO_BINNING_IDX]] }
                | set { fastqs }
         } else {
	        samples = [names, r].transpose()
            channel.from(samples) 
                | map { sample -> [SAMPLE:sample[SAMPLE_IDX],READS:file(sample[READS_IDX])] }
                | set { fastqs }
         }

       emit:
         fastqs
}



workflow _wSRANCBI {
       main:
         MAX_LENGTH=12
         MIN_LENGTH=9
         ACCESSION_ID = 0
         SAMPLE_IDX = 0
         SAMPLE_CONTENT_IDX = 1
         CO_BINNING_IDX = 2

         // Parse TSV file to get access numbers
         accessions = []
         if("sheet" in params.input.SRA.NCBI){
            samples = fetchRunAccessions(params.input.SRA.NCBI.sheet)
            accessions = samples.collect { sample ->  sample[ACCESSION_ID] }
         }

         // check if the number of SRA files is correct and return the correct format 
         FASTQ_LIST = 1 

         // IDs provided via CLI 
         idsFromCLI = []
         if("id" in params.input.SRA.NCBI){
            idsFromCLI = params.input.SRA.NCBI.id.tokenize(" ")
         }

         binGroupFromCLIChannel = Channel.empty()
         if("binGroup" in params.input.SRA.NCBI){
            groupsCLI = params.input.SRA.NCBI.binGroup.tokenize(" ")
            if (idsFromCLI.size() != groupsCLI.size()) {
                error "Mismatch detected: --input.SRA.NCBI.id and --input.SRA.NCBI.binGroup should have the same number of values."
            }
	        coBinningSamplesCLI = [idsFromCLI, groupsCLI].transpose()
            binGroupFromCLIChannel = channel.from(coBinningSamplesCLI)
         }

         Channel.fromSRA(accessions.unique() + idsFromCLI.unique()) | set { foundSRAFiles }

         foundSRAFiles | map { f -> [f[ACCESSION_ID],  Utils.asList(f[FASTQ_LIST])] } 
         | _wCheckSRAFiles 

         Channel.from(accessions) \
		| combine(foundSRAFiles | map { sample -> sample[ACCESSION_ID]} | toList() | toList()) \
  		| filter({ id,idList -> !idList.contains(id) }) | map{ id -> id[ACCESSION_ID] } \
		| set { notFoundAccessions }

        _wCheckSRAFiles.out.passedSamples | map { sample -> [sample.SAMPLE, sample]} 
            | combine(samples, by: SAMPLE_IDX) 
            | map { sample -> sample[SAMPLE_CONTENT_IDX] + ["CO_BINNING":sample[CO_BINNING_IDX]] }   
            | mix(binGroupFromCLIChannel)
            | set { passedSamples }

       emit:
         fastqs = passedSamples
         failedSRAFastqFiles = _wCheckSRAFiles.out.failedSRAIDs
         notFoundAccessions = notFoundAccessions
}


/*
* This method adds the number of samples per co-binned sample, as well as a general flag to indicate whether co-binning should be performed.
*/
def setCoBinning(group_id, samples_list){ 
    // 1. Get the number of elements in the second index
    def total_count = samples_list.size()

    // 2. Add this count to every map in the list
    def updated_list = samples_list.collect { sample_map ->
        return sample_map + [CO_BINNING_COUNT: total_count, DO_CO_BINNING: true]
    } 

    // 3. Return the original structure with the updated list
    return [group_id, updated_list]
}


/*
* This workflow sets the binning specific variables DO_CO_BINNING, CO_BINNING_COUNT and CO_BINNING per sample. 
*/
workflow _wSetCoBinningMetadata {
    take:
        fastqs
    main:
        // Distinguish between samples with and without CO_BINNING column value.
        fastqs | branch { sample ->
            coBinning: sample.containsKey("CO_BINNING")
            singleSample: !sample.containsKey("CO_BINNING")
        } | set { sampleType }

        SAMPLE_IDX=0
        GROUP_IDX=0
        SAMPLE_LIST_IDX=1

        // Set the number of samples per group in every sample array
        // by setting the variables DO_CO_BINNING, CO_BINNING_COUNT and CO_BINNING.
        sampleType.coBinning
            | map { sample -> [sample["CO_BINNING"], sample] }
            | groupTuple(by: GROUP_IDX)
            | map { group_id, samples_list -> setCoBinning(group_id, samples_list) }
            | map { sample ->  sample[SAMPLE_LIST_IDX] } | flatten  | branch { sample -> 
                coBinningSamples: sample["CO_BINNING_COUNT"]>1
                singleBinningSamples: sample["CO_BINNING_COUNT"]==1
            } | set { verifiedSamples }

        // If just a single sample should be processed then indicate that in the metadata.
        sampleType.singleSample  
            | map { sample -> sample + [DO_CO_BINNING: false, CO_BINNING_COUNT:null, CO_BINNING:null] }
            | set { singleSamples }

        // If the CO_BINNING column exists but only one sample is specified, treat that sample as a single sample.
        verifiedSamples.singleBinningSamples
            | map { sample -> sample + [DO_CO_BINNING: false, CO_BINNING_COUNT:null, CO_BINNING:null] }
            | set {verifiedSingleBinningSamples}

        verifiedSamples.coBinningSamples 
            | mix(singleSamples) 
            | mix(verifiedSingleBinningSamples) 
            | set {samples}
    emit:
        samples = samples
}

/*
 *  The input modules defined three input sources: SRA NCBI, generic SRA S3 source that contains a column consisting of SRA IDs 
 *  and a generic source that allows to consume a file containing local path, https or S3 links of paired end or nanopore reads.
 *  SRA NCBI and generic SRA S3 source must contain a column with `ACCESSION` column header and the file for the generic source 
 *  must contain the columns headers `SAMPLE`, `READS1` and `READS2` for paired end and `SAMPLE` and `READS` for nanopore data.
 *  The input workflow allows to process files that are provided via a sample sheet or via CLI.
 * 
 *  In all cases a channel is returned containing values of the format: [TYPE: illumina or ont, SAMPLE:name of the sample, READS1: left read, READS2: right read],
 *  [TYPE: illumina or ont, SAMPLE:name of the sample, READS1: left read, READS2: right read, DO_CO_BINNING: should be co binning done, 
 *  CO_BINNING_COUNT: how many samples should be co-binned, CO_BINNING: name of the group of samples that should be co-binned] 
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
        // Check whether a sample spreadsheet or files are provided on the CLI.
        if("sheet" in params.input.paired) {
            _wSplitReadsSheet().fastqs \
                | map { data -> data["TYPE"] = "ILLUMINA"; data } \
                | set { pairedChannel }
            fastqs | mix(pairedChannel) | set { fastqs }
        }

        if("r1" in params.input.paired || "r2" in params.input.paired) {
            // Make sure that always left and right read is provided and the sample names
            if("r1" !in params.input.paired || "r2" !in params.input.paired || "names" !in params.input.paired){
                error "Missing parameter: --input.paired.r1 and --input.paired.r2 must be provided!"
            }

            _wSplitReadsFiles().fastqs \
              | map { data -> data["TYPE"] = "ILLUMINA"; data } \
              | set { pairedChannel }
            fastqs | mix(pairedChannel) | set { fastqs }
        }
    }

    if("ont" in inputTypes) {
        // Check whether a sample spreadsheet or files are provided on the CLI.
        if("sheet" in params.input.ont) {
            _wOntReadsSheet().fastqs \
                | map { data -> data["TYPE"] = "OXFORD_NANOPORE"; data } \
                | set { ontChannel }

            fastqs | mix(ontChannel) | set { fastqs }
        }

        if("r" in params.input.ont || "names" in params.input.ont) {
          // Make sure that always the reads and sample names are provided 
          if("r" !in params.input.ont || "names" !in params.input.ont){
                error "Missing parameter: --input.ont.r and --input.ont.names must be provided!"
          }
          _wOntReadsFiles().fastqs \
                | map { data -> data["TYPE"] = "OXFORD_NANOPORE"; data } \
                | set { ontChannel }

          fastqs | mix(ontChannel) | set { fastqs }
        }
    }

    // Set co binning specific metadata such as which samples should be co binned.
    fastqs | _wSetCoBinningMetadata
  emit:
    data = _wSetCoBinningMetadata.out.samples 
}
