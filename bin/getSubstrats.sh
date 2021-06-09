zcat $1 | jq '.|.reactions[] | .metabolites | to_entries | .[] | select(.value<0) | .key ' |  tr -d '"'
