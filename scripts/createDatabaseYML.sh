#!/bin/bash

extractedPath="$1"
md5sum="$2"
https="$3"
compressedPath="$4"
s3File="$5"
s3Directory="$6"
s5cmdCommand="$7"
basePath="$8"
yamlToTest="${9}"
yamlKey="${10}"
skipTests="${11}"

if [[ -z "$skipTests" ]] || [[ $skipTests != *"extracted"* ]]; then
database="{extractedDBPath: ${extractedPath} }" \
	yq ${yamlKey} ${yamlToTest} > ${basePath}/extractedPath.yml
fi

if [[ -z "$skipTests" ]] || [[ $skipTests != *"https"* ]]; then
database="{download: {source: ${https}, md5sum: ${md5sum}} }" \
	yq ${yamlKey} ${yamlToTest} > ${basePath}/https.yml
fi

if [[ -z "$skipTests" ]] || [[ $skipTests != *"compressed"* ]]; then
database="{download: {source: ${compressedPath}, md5sum: ${md5sum}} }" \
	yq ${yamlKey} ${yamlToTest} > ${basePath}/compressedPath.yml
fi

if [[ -z "$skipTests" ]] || [[ $skipTests != *"s3File"* ]]; then
database="{download: {source: ${s3File}, md5sum: ${md5sum}, s5cmd: { params: ${s5cmdCommand} }}}" \
	yq ${yamlKey} ${yamlToTest} > ${basePath}/s3File.yml
fi

if [[ -z "$skipTests" ]] || [[ $skipTests != *"s3Directory"* ]]; then
database="{download: {source: ${s3Directory}, md5sum: ${md5sum}, s5cmd: { params: ${s5cmdCommand} }}}" \
	yq ${yamlKey} ${yamlToTest} > ${basePath}/s3Directory.yml
fi

