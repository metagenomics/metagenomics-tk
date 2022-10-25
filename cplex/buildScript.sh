#!/bin/bash

IMAGE=$1
DOCKERFILE_FOLDER="docker folder"

if [[ "$(docker images -q $1 2> /dev/null)" == "" ]]; then
    cd ${DOCKERFILE_FOLDER} && docker build -t ${IMAGE} .
fi
