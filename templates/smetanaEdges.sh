
mkdir results
mkdir failed
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
	tool_out="results/${NUM}_out.tsv";
	failed_out="failed/${NUM}_out.tsv";

	if ! grep -q "n/a" ${smetana_out}; then
		sed "s/all/$V1\t$V2/g" ${smetana_out} | tail -n 1 > ${tool_out}
        else
		sed "s/all/$V1\t$V2/g" ${smetana_out} | tail -n 1 > ${failed_out}
	fi;'

find failed -name '*_out.tsv'  -exec cat {} \; > retry.tsv
find results -name '*_out.tsv'  -exec cat {} \; > computed_out.tsv

# Smetana sometimes fails to compute the MRO or MIP metrics which could also be due to a bad quality of the input MAGs.
# If this is the case it should retry the computation.
mkdir rescue
grep -n "n/a" retry.tsv | while IFS= read -r line ; do  
	NUM=$(echo "$line" | cut -f 1 -d ':'); 
	V1=$(echo "$line" | cut -d$'\t' -f 1 | cut -d$':' -f 2); 
	V2=$(echo "$line" | cut -d$'\t' -f 2);

	# Retry the computation 10 times
	for counter in {1..10}; do 
		echo "Try to rescue $V1 $V2; Retry Counter: $counter ";
		smetana $V1 $V2 !{params.steps.cooccurrence.metabolicAnnotation.additionalParams.smetana} -o rescue/$NUM;
		sed "s/all/$V1\t$V2/g"  rescue/${NUM}_global.tsv | tail -n 1 > rescue/${NUM}_out.tsv ;

		# If we are able to compute the missing value then we can exit the loop
		if ! grep -q "n/a" rescue/${NUM}_global.tsv; then
			echo "Rescued $V1 $V2"; break;
		fi;
	done;
done

# Collect edge results
find rescue -name '*_out.tsv'  -exec cat {} \; > possible_rescued_out.tsv
cat *_out.tsv > edge_attributes.tsv

# Add Header line
sed -i "1 i\V1\tV2\tmedium\tsize\tmip\tmro" edge_attributes.tsv
