def getOutput(SAMPLE, RUNID, TOOL, filename){
    return SAMPLE + '/' + RUNID + '/' + params.modules.fragmentRecruitment.name + '/' + 
         params.modules.fragmentRecruitment.version.major + "." +
         params.modules.fragmentRecruitment.version.minor + "." +
         params.modules.fragmentRecruitment.version.patch +
         '/' + TOOL + '/' + filename
}

process pCovermCount {

    label 'small'

    container "${params.ubuntu_image}"

    publishDir params.output, mode: "${params.publishDirMode}", saveAs: { filename -> getOutput("${sample}", params.runid, "coverm", filename) }

    when:
    run

    input:
      val(covermParams)
      val(run)
      tuple val(sample), file(mapping), file(listOfRepresentatives), val(medianQuality)

    output:
      tuple val("${sample}"), path("${sample}_stats_out/coveredBases.tsv"), emit: mean
      tuple val("${sample}"), path("${sample}_stats_out/metrics.tsv"), emit: metrics
      tuple val("${sample}"), path("foundGenomes.tsv"), optional: true, emit: foundGenomes
      tuple file(".command.sh"), file(".command.out"), file(".command.err"), file(".command.log")

    shell:
    DO_NOT_ESTIMATE_QUALITY = -1 
    MEDIAN_QUALITY=Double.parseDouble(medianQuality)
    percentIdentity = MEDIAN_QUALITY != DO_NOT_ESTIMATE_QUALITY ? \
	" --min-read-percent-identity "+Utils.getMappingIdentityParam(MEDIAN_QUALITY) : " "
    '''
    OUT=!{sample}_stats_out
    mkdir $OUT
    readlink -f !{listOfRepresentatives} > list.txt 
   
    # Create a mapping between file path basename of the file without ending. (/path/test.1.tsv --> test.1) 
    paste -d '\t' list.txt  <(cat list.txt  | rev | cut -d '/' -f 1  | cut -d '.' -f 2- | rev) > mapping.tsv
    
    # Get covered bases
    coverm genome -t !{task.cpus} !{covermParams} -b !{mapping} \
        !{percentIdentity}  \
        --genome-fasta-list list.txt --methods covered_bases --output-file covTmpContent.tsv \

    # Get length
    coverm genome -t !{task.cpus} -b !{mapping} --min-covered-fraction 0  \
        --genome-fasta-list list.txt --methods length --output-file lengthTmpContent.tsv \

    # Join length and covered bases
    join -t$'\t' -1 1 -2 1 covTmpContent.tsv lengthTmpContent.tsv > covLengthTmpContent.tsv

    # Exchange header ad add covered fraction column
    sed -i  -e '1 s/^.*$/SAMPLE\tGENOME\tCOVERED_BASES\tLENGTH/' -e "2,$ s/^/!{sample}\t/g" covLengthTmpContent.tsv  \
                && echo "COVERED_FRACTION" > covLengthTmp.tsv \
		&& awk '(NR>1){ tmp=($3/($4/100)) ; printf"%0.2f\\n", tmp }' covLengthTmpContent.tsv >> covLengthTmp.tsv \
                && paste -d$'\t' covLengthTmpContent.tsv covLengthTmp.tsv > $OUT/coveredBases.tsv || true

    # Run other metrics like RPKM, TPM, ...
    coverm genome  -t !{task.cpus}  -b !{mapping} \
         !{covermParams} !{percentIdentity} \
	--genome-fasta-list list.txt --methods mean trimmed_mean variance length count reads_per_base rpkm tpm \
	| sed -e '1 s/^.*$/SAMPLE\tGENOME\tMEAN\tTRIMMED_MEAN\tVARIANCE\tLENGTH\tREAD_COUNT\tREADS_PER_BASE\tRPKM\tTPM/' \
	| sed -e "2,$ s/^/!{sample}\t/g" > $OUT/metrics.tsv || true

    coveredBasesCutoff=!{params.steps?.fragmentRecruitment?.mashScreen?.coveredBasesCutoff}
    FOUND_GENOMES=foundGenomes.tsv
    for b in $(awk -v coveredBases=${coveredBasesCutoff} '(NR>1){if ($5 > coveredBases) print $2}' $OUT/coveredBases.tsv); do 
        file=$(grep -P "\t$b$" mapping.tsv | cut -f 1);
	echo $(basename $file)  >> ${FOUND_GENOMES} 
    done
    '''

}
