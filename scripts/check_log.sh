LOG=$1

cut -d$'\t' -f 4 ${LOG} | sort | uniq -c | tr -s ' ' | sed 's/^ //g' | grep -P "^2 " | while read -r line ; do
	grep "$(echo $line | cut -d ' ' -f 2-)"  ${LOG} | grep -q COMPLETED || (echo "$?"; exit 1);
done