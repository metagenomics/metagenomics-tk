GENOMES=$(paste -d$'\t' <(ls -1 genome1*) <(ls -1 genome2*));
while IFS= read -r genomes_to_compare; do 
    GENOME1=$(echo "$genomes_to_compare" | cut -d$'\t' -f 1)
    GENOME2=$(echo "$genomes_to_compare" | cut -d$'\t' -f 2)

    dir=${GENOME1}_vs_${GENOME2}_dir
    mkdir $dir
    cd $dir
    ln -s ../$GENOME1 ${GENOME1}.fa
    ln -s ../$GENOME2 ${GENOME2}.fa
    average_nucleotide_identity.py !{params.dereplication.pyani_parameters} --force -i . -o out &> error.log
    ANI1=$(cat out/*_percentage_identity.tab | cut -f 2 | tail -n 1) 
    ANI2=$(cat out/*_percentage_identity.tab | cut -f 3 | head -n 2 | tail -n 1)
    ANI=$(echo | awk -v ani1="$ANI1" -v ani2="$ANI2" '{  print (ani1 + ani2) / 2 }')

    GENOME1_PATH=$(readlink -f ../$GENOME1)
    GENOME2_PATH=$(readlink -f ../$GENOME2)

    ANI_RESULT_UPPER="$GENOME1_PATH\t$GENOME2_PATH\t$ANI2"
    echo -e $ANI_RESULT_UPPER > ${GENOME1}_vs_${GENOME2}.tsv.upper

    ANI_RESULT_LOWER="$GENOME1_PATH\t$GENOME2_PATH\t$ANI1"
    echo -e $ANI_RESULT_LOWER > ${GENOME1}_vs_${GENOME2}.tsv.lower

    ANI_RESULT="$GENOME1_PATH\t$GENOME2_PATH\t$ANI"
    echo -e $ANI_RESULT > ${GENOME1}_vs_${GENOME2}.tsv
    cd -

done <<< $GENOMES
DIRECTORY=$(mktemp -d --suffix=.out  -p .)
cat *_dir/out/*_percentage_identity.tab > ${DIRECTORY}/matrix.tab
cat *_dir/*.tsv > ${DIRECTORY}/out.tsv 
cat *_dir/*.tsv.upper > ${DIRECTORY}/out.tsv.upper
cat *_dir/*.tsv.lower > ${DIRECTORY}/out.tsv.lower
