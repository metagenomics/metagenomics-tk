include { pProdigal ; pHmmSearch } from '../annotation/module'

/**
 * MAGScoT - Run MAGScoT binning refinement
 * @param sample: Sample name
 * @param contigMaps: Individual contig to bin files from binning algorithms
 * @param allHits: GTDBtk single copy marker identification of all contigs
 * @param contigs: Path to contigs
 * @return: Files containing scores, bin-contig mappings, new bins, and not binned contigs
 *
 * This process runs the MAGScoT R script on the input contig maps
 * to compare and refine the binning of multiple other binners.
 * It outputs several files containing scores, bin-contig mappings, bins, and not binned contigs.
 **/
process pMAGScoT {

    container "${params.magscot_image}"

    tag "Sample: ${sample}"

    label 'small'

    executor "local"

    publishDir params.output, mode: "${params.publishDirMode}", saveAs: { filename ->
        Output.getOutput("${sample}", params.runid, "refinement/magscot", params.modules.binning, filename)
    }

    // Override default MAGScoT container entrypoint, so that the "RScript magscot" call that is normally run
    // does not clash with the Nextflow process call "/bin/bash"
    containerOptions params.apptainer ? "" : ' --entrypoint "" '

    when params.steps.containsKey("binRefinement") && params.steps.binRefinement.containsKey("magscot")

    input:
    tuple val(sample), file(contigMaps), file(allHits), path(contigs)

    output:
    tuple val("${sample}"), file("${sample}_MagScoT.*"), optional: true, emit: scores
    tuple val("${sample}"), file("${sample}_bin_contig_mapping.tsv"), optional: true, emit: binContigMapping
    tuple val("${sample}"), path("${sample}_bin.*.fa", arity: '0..*'), emit: bins
    tuple val("${sample}"), file("${sample}_notBinned.fa"), optional: true, emit: notBinned
    tuple file(".command.sh"), file(".command.out"), file(".command.err"), file(".command.log")

    shell:
    '''
    # Once in a blue moon Nextflow leaves the header in the file
    # Failsafe to remove header from contigMaps file if it exists
    sed -i '/^BIN_ID\tCONTIG\tBINNER$/d' !{contigMaps}
    Rscript /opt/MAGScoT.R !{params.steps?.binRefinement?.magscot?.additionalParams} -i !{contigMaps} --hmm !{allHits} -o !{sample}_MagScoT

    # Create a new binning file according to the naming convention
    echo "Converting MAGScoT binning according to the naming convention"
    echo -e 'BIN_ID\tCONTIG\tBINNER' > !{sample}_bin_contig_mapping.tsv

    # Use head to get the first line, awk to print the first field (the filename),
    # and sed to remove everything up to and including the last dot (.) in the filename, leaving only the file extension.
    # Export the variable to use it in the xargs commands subshell further down
    export EXT=$(head -n 1 !{contigMaps} | awk '{print $1}' | sed 's/.*\\.//')

    # Remove the header and pipe the remaining lines to xargs to run the script line by line
    sed 1d !{sample}_MagScoT.refined.contig_to_bin.out | xargs -n 2 sh -c '
        # Get the first column and separate the number, remove the leading zeros
        binID=$(echo $0 | rev | cut -d"_" -f1 | rev | sed 's/^0*//')
        CONTIG=$1
        # Create a file with the contigs for each bin for reconstruction with seqkit
        echo $CONTIG >> !{sample}_bin.$binID.lst
        echo "!{sample}_bin.$binID.$EXT\t$CONTIG\tMAGScot" >> !{sample}_bin_contig_mapping.tsv
    '
    echo "Done converting"
    # Reconstructing the bins from the converted MAGScoT output list
    echo "Reconstructing bins"
    for bin in !{sample}_bin.*.lst; do
        seqkit grep -f $bin !{contigs} -o ${bin%.lst}.fa.tmp
    done
    echo "Done reconstructing bins"

    # Rename the bins to the naming convention
    for bin in $(find * -name "*bin*.fa.tmp"); do
    	BIN_NAME="$(basename ${bin})"
    	echo "Renaming bin: ${BIN_NAME}"

    	# Get id of the bin (e.g get 2 of the bin SAMPLEID_bin.2.fa.tmp)
    	ID=$(echo ${BIN_NAME} | rev | cut -d '.' -f 3 | rev)
    	echo "Bin id: ${ID}"

    	# Append bin id to every header
    	seqkit replace  -p '(.*)' -r "\\${1} MAG=${ID}" $bin > ${BIN_NAME%.tmp}
    	# Remove temporary file
    	rm $bin
    done

    # Creating un-binned contigs file
    seqkit grep -f <(cat !{sample}_bin.*.lst) -v !{contigs} -o !{sample}_notBinned.fa.tmp
    seqkit replace -p '(.*)' -r "\\${1} MAG=NotBinned" !{sample}_notBinned.fa.tmp > !{sample}_notBinned.fa
    '''
}


process pBinette {

    container "${params.binette_image}"

    tag "Sample: ${sample}"

    label 'small'

    containerOptions Utils.getDockerMount(params.steps?.binRefinement?.binette?.database, params, apptainer=params.apptainer) + (params.apptainer ? "" : Utils.getDockerNetwork()) + "-u \$(id -u):\$(id -g)" 

    publishDir params.output, mode: "${params.publishDirMode}", saveAs: { filename ->
        Output.getOutput("${sample}", params.runid, "refinement/binette", params.modules.binning, filename)
    }

    input:
    tuple val(sample), path(contigMaps, name: "contigMaps/contigMap*.tsv"), path(contigs)

    output:
    tuple val("${sample}"), file("${sample}_bin_contig_mapping.tsv"), optional: true, emit: binContigMapping
    tuple val("${sample}"), path("${sample}_bin.*.fa", arity: '0..*'), emit: bins
    tuple val("${sample}"), file("${sample}_notBinned.fa"), optional: true, emit: notBinned
    tuple file(".command.sh"), file(".command.out"), file(".command.err"), file(".command.log")

    script:
    S5CMD_PARAMS=params?.steps?.binRefinement?.binette?.database?.download?.s5cmd?.params ?: "" 
    DOWNLOAD_LINK=params?.steps?.binRefinement?.binette?.database?.download?.source ?: ""
    MD5SUM=params.steps?.binRefinement?.binette?.database?.download?.md5sum ?: ""
    EXTRACTED_DB=params.steps?.binRefinement?.binette?.database?.extractedDBPath ?: ""
    S3_checkm2_ACCESS=params?.steps?.binRefinement?.binette?.database?.download?.s5cmd && S5CMD_PARAMS.indexOf("--no-sign-request") == -1 ? "\$S3_checkm2_ACCESS" : ""
    S3_checkm2_SECRET=params?.steps?.binRefinement?.binette?.database?.download?.s5cmd && S5CMD_PARAMS.indexOf("--no-sign-request") == -1 ? "\$S3_checkm2_SECRET" : ""
    """
    if [ -z "${EXTRACTED_DB}" ]
    then
      DATABASE=${params.polished.databases}/checkm2
      LOCK_FILE=\${DATABASE}/lock.txt

      if [ ! -z "${S3_checkm2_ACCESS}" ]
      then
        export AWS_ACCESS_KEY_ID=${S3_checkm2_ACCESS}
        export AWS_SECRET_ACCESS_KEY=${S3_checkm2_SECRET}
      fi

      # Download checkm database if necessary
      mkdir -p \${DATABASE}
      flock \${LOCK_FILE} concurrentDownload.sh --output=\${DATABASE} \
        --link=${DOWNLOAD_LINK} \
        --httpsCommand="wgetStatic --no-check-certificate -O checkm2.tar.gz ${DOWNLOAD_LINK} && tar -xzvf checkm2.tar.gz && rm checkm2.tar.gz" \
        --s3FileCommand="s5cmd ${S5CMD_PARAMS} cp --concurrency ${task.cpus} ${DOWNLOAD_LINK} checkm2.tar.gz && tar -xzvf checkm2.tar.gz && rm checkm2.tar.gz" \
        --s3DirectoryCommand="s5cmd ${S5CMD_PARAMS} cp --concurrency ${task.cpus} ${DOWNLOAD_LINK} . " \
        --s5cmdAdditionalParams="${S5CMD_PARAMS}" \
        --localCommand="tar -xzvf ${DOWNLOAD_LINK}" \
        --expectedMD5SUM=${MD5SUM}
     
      export CHECKM2DB=\$(find \${DATABASE} -name "*.dmnd")
    else
      export CHECKM2DB=\$(find ${EXTRACTED_DB} -name "*.dmnd")
    fi

    mkdir binetteInput
    for map in \$(find contigMaps ! -type d); do 
      csvtk cut -t -f CONTIG,BIN_ID \${map} | tail -n +2 > binetteInput/\$(basename \${map})
    done

    binette ${params.steps.binRefinement.binette.additionalParams} --contig2bin_tables binetteInput/* --contigs ${contigs} --checkm2_db \${CHECKM2DB} 

    BIN_CONTIG_MAPPING=${sample}_bin_contig_mapping.tsv
    echo -e "BIN_ID\tCONTIG\tBINNER" > \${BIN_CONTIG_MAPPING}
    for bin in \$(find results/final_bins -name "bin*.fa"); do
      BIN_NUM=\$(echo "\$bin" | grep -oE 'bin[0-9]+' | grep -oE '[0-9]+')

      BIN_NAME="${sample}_bin.\${BIN_NUM}.fa"

      # Get id of the bin (e.g get 2 of the bin SAMPLEID_bin.2.fa)
      ID=\$(echo \${BIN_NAME} | rev | cut -d '.' -f 2 | rev)

      # Append bin id to every header
      seqkit replace  -p '(.*)' -r "\\\${1} MAG=\${ID}" \$bin > \${BIN_NAME}

      # Create bin to contig mapping and add the used binner to each line
      grep ">" \${bin} | sed 's/>//g' \
        | sed "s/^/\${BIN_NAME}\t/g;s/\$/\tmetabat/" >> \${BIN_CONTIG_MAPPING}
    done

    # return not binned fasta files
    BINNED_IDS=binned.tsv
    NOT_BINNED=${sample}_notBinned.fa
    grep -h ">" \$(find results/final_bins -name "bin*.fa") | tr -d ">" > \${BINNED_IDS}
    if [ -s \${BINNED_IDS} ]; then
      # Get all not binned Ids
      seqkit grep -vf \${BINNED_IDS} ${contigs} \
        | seqkit replace  -p '(.*)' -r "\\\${1} MAG=NotBinned" > \${NOT_BINNED}
    else
      seqkit replace  -p '(.*)' -r "\\\${1} MAG=NotBinned" ${contigs} > \${NOT_BINNED}
    fi
    """
}



process pAsTools {

    container "${params.agtools_image}"

    tag "Sample: ${sample}"

    label 'small'

    containerOptions params.apptainer ? "" : Utils.getDockerNetwork() + "-u \$(id -u):\$(id -g)" 

    publishDir params.output, mode: "${params.publishDirMode}", saveAs: { filename ->
        Output.getOutput("${sample}", params.runid, "refinement/binSpreader", params.modules.binning, filename)
    }

    input:
    tuple val(sample), path(fastg), val(kmerSize)

    output:
    tuple val("${sample}"), file("${sample}_assembly.gfa"), optional: true, emit: gfa 
    tuple file(".command.sh"), file(".command.out"), file(".command.err"), file(".command.log")

    script:
    """
    agtools fastg2gfa --graph  ${fastg} --ksize ${kmerSize} --output  ${sample}_assembly.gfa
    """
}

process pBinSpreader {

    container "${params.metaspades_image}"

    tag "Sample: ${sample}"

    label 'small'

    containerOptions params.apptainer ? "" : Utils.getDockerNetwork()

    publishDir params.output, mode: "${params.publishDirMode}", saveAs: { filename ->
        Output.getOutput("${sample}", params.runid, "refinement/binSpreader", params.modules.binning, filename)
    }

    input:
    tuple val(sample), path(contigMaps), path(contigs), path(gfa), val(maxKmer), path(paths), path(headerMapping)

    output:
    tuple val("${sample}"), file("${sample}_bin_contig_mapping.tsv"), optional: true, emit: binContigMapping
    tuple val("${sample}"), path("bins/${sample}_bin.*.fa", arity: '0..*'), emit: bins
    tuple val("${sample}"), file("bins/${sample}_notBinned.fa"), optional: true, emit: notBinned
    tuple file(".command.sh"), file(".command.out"), file(".command.err"), file(".command.log")

    script:
    """
    csvtk replace -t -f CONTIG -p "^(.+)\$" -r '{kv}' -k <(awk '{print \$2 "\t" \$1}' ${headerMapping}) ${contigMaps} \
        | csvtk cut -t -f CONTIG,BIN_ID \
        | tail -n +2 > binSpreaderInputMap.tsv

    binspreader ${gfa} binSpreaderInputMap.tsv out -t ${task.cpus} --paths ${paths} ${params.steps.binRefinement.binSpreader.additionalParams}

    csvtk replace -H -t -f 1 -p "^(.+)\$" -r '{kv}' -k ${headerMapping} out/binning.tsv  \
       | csvtk cut -t -f 2,1 > renamed_contig_binning.tsv

    OUTDIR="bins"

    mkdir -p "\$OUTDIR"

    cut -f1 renamed_contig_binning.tsv | sort -u | while read -r bin; do

        # extract IDs belonging to this bin
        awk -F'\t' -v b="\$bin" '\$1==b {print \$2}' renamed_contig_binning.tsv > "\$OUTDIR/\${bin}.ids.txt"

        # fetch those sequences from the contigs file
        seqkit grep -f "\$OUTDIR/\${bin}.ids.txt" ${contigs} > "\$OUTDIR/\${bin}"

         rm "\$OUTDIR/\${bin}.ids.txt"
    done

    mkdir -p mapping

    BIN_CONTIG_MAPPING=mapping/${sample}_bin_contig_mapping.tsv
    sed '1i BIN_ID\tCONTIG\tBINNER' renamed_contig_binning.tsv  \
     | sed '2,\$s/\$/\tBINSPREADER/' > \${BIN_CONTIG_MAPPING}

    cut -f1 renamed_contig_binning.tsv | sort -u > "\$OUTDIR/binned_ids.txt"

    seqkit grep -v -f "\$OUTDIR/binned_ids.txt" ${contigs} > "\$OUTDIR/${sample}_notBinned.fa"
    rm "\$OUTDIR/binned_ids.txt"  
    """
}


workflow _wMAGScoT {
    take:
    contigs
    binContigMapping

    main:
        SAMPLE_IDX = 0
        CONTIG_MAPPING_IDX = 1
        pProdigal(contigs)
        pHmmSearch(pProdigal.out.prodigal_faa)
        binContigMapping
            | collectFile(keepHeader: false) { item -> ["${item[SAMPLE_IDX]}", item[CONTIG_MAPPING_IDX].text] }
            | map { f -> [file(f).name, f] }
            | join(pHmmSearch.out.allhits, by: SAMPLE_IDX)
            | join(contigs, by: SAMPLE_IDX)
            | set { magscot_input }

        pMAGScoT(magscot_input)
        pMAGScoT.out.bins | set { bins }
        pMAGScoT.out.notBinned | set { notBinned }
        pMAGScoT.out.binContigMapping | set { binContigMapping }
    emit:
    bins = bins
    notBinned = notBinned
    binContigMapping = binContigMapping
}

workflow wRefinementList {

    take:
    contigs
    binContigMapping
    fastg
    gfa
    paths
    headerMapping

    main:
    _wRefinement(contigs, binContigMapping, fastg, gfa, paths, headerMapping)

    emit:
    bins = _wRefinement.out.bins
    notBinned = _wRefinement.out.notBinned
    binContigMapping = _wRefinement.out.binContigMapping
}


/*
*
* This workflow takes an input_reads channel as input with the following format [SAMPLE, READS PAIRED, READS UNPAIRED]
* and a contigs channel with the format [SAMPLE, CONTIGS]
*
*/
workflow _wRefinement {
    take:
    contigs
    binContigMapping
    fastg
    gfa
    paths
    headerMapping

    main:
    SAMPLE_IDX = 0

    bins = channel.empty()
    notBinned = channel.empty()
    // Only use MAGScoT bins if the user has selected the refinement step
    if (params.steps.containsKey("binRefinement") && params.steps.binRefinement.containsKey("magscot")) {
        _wMAGScoT(contigs, binContigMapping)
        _wMAGScoT.out.bins | set { bins }
        _wMAGScoT.out.notBinned | set { notBinned }
        _wMAGScoT.out.binContigMapping | set { binContigMapping }
    } else if (params.steps.containsKey("binRefinement") && params.steps.binRefinement.containsKey("binette")){
        SAMPLE_IDX = 0
        CONTIG_MAPPING_IDX = 1
        binContigMapping
            | groupTuple(by: SAMPLE_IDX)
            | join(contigs, by: SAMPLE_IDX)
            | pBinette

        pBinette.out.bins | set { bins }
        pBinette.out.notBinned | set { notBinned }
        pBinette.out.binContigMapping | set { binContigMapping }
    }

    if (params.steps.containsKey("binRefinement") && params.steps.binRefinement.containsKey("binSpreader")){
        binContigMapping
         | combine(contigs, by: SAMPLE_IDX)
         | combine(gfa, by: SAMPLE_IDX)
         | combine(paths, by: SAMPLE_IDX)
         | combine(headerMapping, by: SAMPLE_IDX)
         | pBinSpreader

        pBinSpreader.out.bins | set { bins }
        pBinSpreader.out.notBinned | set { notBinned }
        pBinSpreader.out.binContigMapping | set { binContigMapping }
    }

    emit:
    bins = bins
    notBinned = notBinned
    binContigMapping = binContigMapping
}
