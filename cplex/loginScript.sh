#!/bin/bash

IMAGE=$1

if [[ "$(docker images -q $1 2> /dev/null)" == "" ]]; then
        echo "your password" | docker login --username pbelmann --password-stdin  && docker pull ${IMAGE}
fi

