# Create Bash array
sraids=($(echo "!{sraids}" | sed -e 's/\[//g' -e 's/]//g' -e 's/,//g'))

# Set header in output files
FOUND_FILE="found.csv"
NOT_FOUND_FILE="not_found.csv"

echo "study_id,run_accession,instrument" > ${FOUND_FILE}
echo "run_accession" > ${NOT_FOUND_FILE}

for sraid in "${sraids[@]}"; do
    entry=""

    # skip retrieval from DB 
    if [[ -z "!{skipDB}" ]]; then
        entry=$(zstd -dc !{projectDir}/assets/modules/input/sra/sra_db.csv.zst | grep -e ",${sraid}," -e "^${sraid}," | cat)
    fi

    # Differentiate between found and not found DB entry
    if [[ -z "$entry" ]]; then
          echo ${sraid} >> ${NOT_FOUND_FILE}
    else
          echo ${entry} >> ${FOUND_FILE}
    fi
done
