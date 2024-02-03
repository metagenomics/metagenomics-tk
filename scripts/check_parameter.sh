OUTPUT=$1

EXPECTED_MEMORY=2000
MEMORY=$(find  $OUTPUT -name 'params_*'  | head -n 1 | xargs -I {}  yq '.resources | .highmemLarge | .memory' {})

if [ "$MEMORY" -eq "$EXPECTED_MEMORY" ]; then
	exit 0
else
	exit 1
fi
