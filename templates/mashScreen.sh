set -o pipefail

mash screen -p !{task.cpus} !{params?.steps?.fragmentRecruitment?.mashScreen?.additionalParams.mashScreen} !{sketch} !{pairedReads} !{singleReads} \
	| sed "s/^/!{sample}\t/g" >> mash_screen.tsv

set +o pipefail

awk -v dist=!{params?.steps?.fragmentRecruitment?.mashScreen?.mashDistCutoff} \
    -v hash=!{params?.steps?.fragmentRecruitment?.mashScreen?.mashHashCutoff} \
    '{ if ($2 >= dist && $3 >= hash) print $7 }' <(cat mash_screen.tsv | sed "s/\//\t/g") > selected_genomes.tsv
