#!/bin/bash

# This script creates multiple configuration files for testing database input configurations and
# runs nextflow using them as input.

# Example: /vol/spool/peter/plsdb
extractedPath="$1"
# Example: caa846fbb689deba0e35ef559793b521
md5sum="$2"
# Example: https://openstack.cebitec.uni-bielefeld.de:8080/databases/plsdb.zip
https="$3"
# Example: /vol/spool/peter/plsdb.zip
compressedPath="$4"
# Example: s3://databases/plsdb.zip
s3File="$5"
# Example: "s3://databases/plsdb/"
s3Directory="$6"
# Example: '--retry-count 30 --no-verify-ssl --endpoint-url https://openstack.cebitec.uni-bielefeld.de:8080'
s5cmdCommand="$7"
# Example: /vol/path/to/credentials
s5cmdKey="$8"
# Example: "/vol/spool/peter/meta-omics-toolkit/generated_yamls/"
basePath="$9"
# Example: example_params/plasmid.yml
yamlToTest="${10}"
# Example: ".steps.plasmid.PLSDB.database=env(database)"
yamlKey="${11}"
# Example: "./scripts/test_plasmids.sh"
scriptToTest="${12}"
# Example: "/vol/spool/peter/meta-omics-toolkit/plasmid_yaml_database_tests"
workDirBasePath="${13}"
# Example: https,s3File
skipTests="${14}"
# Example: true
deleteScratchDir="${15:-yes}"

mkdir -p ${basePath}

# Create all possible yaml configuration files
bash ./scripts/createDatabaseYML.sh "${extractedPath}" "${md5sum}" "${https}" \
	"${compressedPath}" "${s3File}" "${s3Directory}"  \
	"${s5cmdCommand}" "${s5cmdKey}" \
	"${basePath}" "${yamlToTest}" "${yamlKey}" "${skipTests}"

# Run tests on all configuration files
for yml in ${basePath}/* ; do 
	echo "Configuration File ${yml}"; 
        workDir=${workDirBasePath}/work_$(basename ${yml})
	mkdir -p ${workDir}
	databaseDir="/vol/scratch/databases_$(cat /proc/sys/kernel/random/uuid)/"
	bash ${scriptToTest} " --databases=${databaseDir} " ${yml} "${workDir}" "slurm" || exit 1
	# Cleanup scratch directory
        if [[ "$deleteScratchDir" == "yes" ]]; then
		sinfo -o " %n %t"   | tail -n +2 | sed 's/^ //g' | cut -d ' ' -f 1 | grep -v "dummy" | xargs -P 10 -I {}  ssh -tt ubuntu@{} " sudo rm -rf ${databaseDir} "
	fi
done
