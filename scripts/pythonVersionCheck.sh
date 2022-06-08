#!/bin/bash

version=$(python3 -V 2>&1 | grep -Po '(?<=Python )(.+)')
if [[ -z "$version" ]]
then
	           echo "No Python3 installation found!" 
fi

MINIMUM_PYTHON_VERSION=3.9.0
MINIMUM_PYTHON_VERSION_PARSED=$(echo ${MINIMUM_PYTHON_VERSION} | tr -d '.')

parsedVersion=$(echo "${version//./}")
if [[ "$parsedVersion" -ge "${MINIMUM_PYTHON_VERSION_PARSED}" ]]
then
	           exit 0
	   else
                   echo "Minimum Python version required: ${MINIMUM_PYTHON_VERSION}" 
                   exit 1
fi

