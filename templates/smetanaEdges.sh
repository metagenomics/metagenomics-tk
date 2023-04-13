
nl  -nln !{edges} \
	| cut -f 1 \
	| xargs -i -P 2 bash -c '
	NUM_INPUT="{}"; 
	NUM=$(echo -n ${NUM_INPUT});

	V=$(sed "$NUM!d" !{edges});
	V1=$(echo $V | cut -f 1 -d " ");
	V2=$(echo $V | cut -f 2 -d " "); 

	# Run SMETANA
	smetana $V1 $V2 !{params.steps.cooccurrence.metabolicAnnotation.additionalParams.smetana} -o $NUM;

	smetana_out="${NUM}_global.tsv";
	tool_out="${NUM}_out.tsv";

	sed "s/all/$V1\t$V2/g" ${smetana_out} | tail -n 1 > ${tool_out} '

# Collect edge results
cat *_out.tsv > edge_attributes.tsv

# SMETANA sometimes fails to compute the MRO. If this is the case it should report
# the exit status and retry the computation
if grep -q "n/a" edge_attributes.tsv; then
        echo "There is at least one MRO set n/a in the output of Smetana."
	exit 1
fi

# Add Header line
sed -i "1 i\V1\tV2\tmedium\tsize\tmip\tmro" edge_attributes.tsv
