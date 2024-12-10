
BINS_DIR=bins

binPrefix=$(find $BINS_DIR -name "*_bin.*.fa" -exec readlink -f {} \; \
				            | tail -n 1 \
				            | rev | cut -f 1 -d '/' \
			                    | rev | cut -d '.' -f 1 \
					    | sed 's/$/./g')

annotatedbins2json -checkm-tsv  !{checkm} \
	-gtdbtk-tsvs !{gtdbtk} \
	-json-gz !{sample}.bins.json.gz \
        -bin-id-prefix ${binPrefix}
