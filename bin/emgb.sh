#!/bin/bash
set -e

VERSION=0.4.0

while [ $# -gt 0 ]; do
	  case "$1" in
	    --output=*) OUTPUT_PATH="${1#*=}"
	    ;;
	    --runid=*) RUN_ID="${1#*=}"
	    ;;
	    --binsdir=*) BINS_DIR="${1#*=}"
	    ;;
	    --name=*) NAME="${1#*=}"
	    ;;
	    --db=*) DB="${1#*=}"
	    ;;
	    --workdir=*) WORK_DIR="${1#*=}"
	    ;;
            --type=*) TYPE="${1#*=}"
	    ;;
            --titles=*) TITLES="${1#*=}"
	    ;;
	    --blastdb=*) BLAST_DB="${1#*=}"
	    ;;
	    --version) VERSION_CHECK=1
	    ;;
	    --debug) DEBUG_CHECK=1
	    ;;
	    --help) HELP_CHECK=1
	    ;;
	    *)
              printf "***************************\n"
              printf "* Error: Invalid argument.*\n"
	      printf "***************************\n"
	      exit 1
	    ;;
	  esac
          shift
done



function getGenes {
	nr=$(find $OUTPUT_PATH/$RUN_ID/annotation/ -name "*.${BLAST_DB}.blast.tsv" -exec readlink -f {} \;  | sed 's/^/ --blast-tab /g')
	tax=$(find $OUTPUT_PATH/$RUN_ID/annotation/ -name "*.taxonomy.tsv" -exec readlink -f {} \; | sed 's/^/ -mmseqs-lineage /g')
	ffn=$(find $OUTPUT_PATH/$RUN_ID/annotation -name "*.ffn.gz" -exec readlink -f {} \; | sed 's/^/ -ffn /g')
	gff=$(find $OUTPUT_PATH/$RUN_ID/annotation -name "*.gff.gz" -exec readlink -f {} \; | sed 's/^/ -gff /g')
	faa=$(find $OUTPUT_PATH/$RUN_ID/annotation -name "*.faa.gz" -exec readlink -f {} \; | sed 's/^/ -faa /g')
	kegg=$(find $OUTPUT_PATH/$RUN_ID/annotation/ -name "*.kegg.blast.tsv" -exec readlink -f {} \; | sed 's/^/ -kegg-blast-tab /g')
	db=$DB
	json=" -json-gz $(pwd)/${NAME}.genes.json.gz "
	name=" -dataset-name ${NAME} "

	titles_mount=""
	titles=""
	if [ ! -z "$TITLES" ]
	then
                titles=" -title-tsv $TITLES "
		titles_mount=" -v $TITLES:$TITLES "
	fi

	cmd="$nr $kegg $gff $faa $ffn $tax $bins $db $json $name $titles "

	if [ ! -z "$DEBUG_CHECK" ]
	then
		echo $cmd
	fi

	docker run -i $DBMOUNT -v $(pwd):$(pwd) ${titles_mount}  -v $WORK_DIR:$WORK_DIR -v ${OUTPUT_PATH}:${OUTPUT_PATH} quay.io/emgb/annotatedgenes2json:2.5.0 $cmd
}



function getContigs {
	contigs=$(find  $OUTPUT_PATH/$RUN_ID/assembly${TYPE}/ -name "*_contigs.fa.gz" -exec readlink -f {} \; | sed 's/^/ -fasta /g')
        name=" -sample-names ${NAME} "	
	bam=$(find  $OUTPUT_PATH/$RUN_ID/binning${TYPE}/ -name "*.bam" -exec readlink -f {} \; | sed 's/^/  -sample-bam-files  /g')
	json=" -json-gz $(pwd)/${NAME}.contigs.json.gz "


	cmd="$contigs $name $bam $bins $json"

	if [ ! -z "$DEBUG_CHECK" ]
	then
		echo $cmd
	fi

	docker run -i -v $(pwd):$(pwd) -v $WORK_DIR:$WORK_DIR -v ${OUTPUT_PATH}:${OUTPUT_PATH} quay.io/emgb/annotatedcontigs2json:2.2.2 $cmd
}


function getBins {
	checkm=$(find $OUTPUT_PATH/$RUN_ID/magAttributes/*/checkm*/ -name "*_checkm*_*.tsv" -exec readlink -f {} \; | sed 's/^/ -checkm-tsv /g')
	gtdbtk=$(find $OUTPUT_PATH/$RUN_ID/magAttributes/*/gtdb/ -name "*.summary.tsv" -exec readlink -f {} \; | sed 's/^/ -gtdbtk-tsvs /g')
	bins=$(find $BINS_DIR -name "*_bin.*.fa" -exec readlink -f {} \; | tail -n 1 | rev | cut -f 1 -d '/' | rev | cut -d '.' -f 1 | sed 's/^/  -bin-id-prefix /g')
	json=" -json-gz $(pwd)/${NAME}.bins.json.gz "
 
	cmd="$checkm $gtdbtk $bins. $json"

	if [ ! -z "$DEBUG_CHECK" ]
	then
		echo $cmd
	fi

	docker run -i $DBMOUNT -v $(pwd):$(pwd) -v $WORK_DIR:$WORK_DIR -v ${OUTPUT_PATH}:${OUTPUT_PATH} quay.io/emgb/annotatedbins2json:2.3.0 $cmd
}


help()
{
	echo "  --output    -- toolkit specific sample output directory"
	echo "              -- (e.g. fullPipelineOutput/SAMPLE) "
	echo "  --runid     -- run id of the specific toolkit output"
	echo "              -- (e.g. X in the following example path fullPipelineOutput/SAMPLE/X/binning/)                       "
	echo "  --binsdir   -- directory of bins. If bin refinement was executed then the bin refinement output should be used."
	echo "              -- (e.g. --binsdir=fullPipelineOutput/DRR066656/1/binning/0.4.0/metabat)"
	echo "  --blastdb   -- Blast output that should be exported to emgb"
	echo "              -- (e.g. the folder name of BLAST_DB: output/test1/1/annotation/0.3.0/mmseqs2/BLAST_DB)"
	echo "              -- (Examples: bacmet20_predicted, ncbi_nr)"
	echo "  --db        -- emgb specific kegg database"
	echo "  --name      -- sample name, e.g. the SAMPLE in the paths above"
	echo "  --titles    -- Custom .tsv file mapping IDs to Subject Titles extracted from FASTA (optional)."
	echo "              -- The following link explains how to create the corresponding tsv file."
	echo "              -- https://gitlab.ub.uni-bielefeld.de/cmg/emgb/annotatedgenes2json#use-cases"
	echo "  --type      -- if other then Illumina: ONT/Hybrid"
	echo "  --workdir   -- absolute path to Nextflow work directory"
	echo "  --help      -- help page"
	echo "  --debug     -- print commands before running"
	echo "  --version   -- version of this script"
}

if [ ! -z "$VERSION_CHECK" ]
then
	echo "$VERSION"
	exit 0
fi

if [ ! -z "$HELP_CHECK" ]
then
	help
	exit 0
fi

OUTPUT_PATH=$(readlink -f $OUTPUT_PATH)
bins=" -bins-dir $(readlink -f $BINS_DIR)"

if [ -z "$DB" ]
then
	DBMOUNT=""
   	DB=" -ci "
else
	DBMOUNT=" -v $DB:$DB "
	DB=" -db ${DB} "
fi

getGenes
getContigs
getBins
