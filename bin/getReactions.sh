cat $1 | jq '.|.reactions[] |.name ' | sort
