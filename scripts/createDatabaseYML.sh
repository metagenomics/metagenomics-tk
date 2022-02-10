#!/bin/bash

extractedPath="$1"
md5sum="$2"
https="$3"
compressedPath="$4"
s3File="$5"
s3Directory="$6"
s5cmdCommand="$7"
s5cmdKey="$8"
basePath="$9"
yamlToTest="${10}"
yamlKey="${11}"
skipTests="${12}"

if [[ "extracted" != *"$skipTests"* ]]; then
database="{extractedDBPath: ${extractedPath} }" \
	yq ${yamlKey} ${yamlToTest} > ${basePath}/extractedPath.yml
fi

if [[ "https" != *"$skipTests"* ]]; then
database="{download: {source: ${https}, md5sum: ${md5sum}} }" \
	yq ${yamlKey} ${yamlToTest} > ${basePath}/https.yml
fi

if [[ "compressed" != *"$skipTests"* ]]; then
database="{download: {source: ${compressedPath}, md5sum: ${md5sum}} }" \
	yq ${yamlKey} ${yamlToTest} > ${basePath}/compressedPath.yml
fi

if [[ "s3File" != *"$skipTests"* ]]; then
database="{download: {source: ${s3File}, md5sum: ${md5sum}, s5cmd: { params: ${s5cmdCommand}, keyfile: ${s5cmdKey} }}}" \
	yq ${yamlKey} ${yamlToTest} > ${basePath}/s3File.yml
fi

if [[ "s3Directory" != *"$skipTests"* ]]; then
database="{download: {source: ${s3Directory}, md5sum: ${md5sum}, s5cmd: { params: ${s5cmdCommand}, keyfile: ${s5cmdKey} }}}" \
	yq ${yamlKey} ${yamlToTest} > ${basePath}/s3Directory.yml
fi




