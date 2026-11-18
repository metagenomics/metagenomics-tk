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
# Example: "/vol/spool/peter/meta-omics-toolkit/generated_yamls/"
basePath="$8"
# Example: example_params/plasmid.yml
yamlToTest="${9}"
# Example: ".steps.plasmid.PLSDB.database=env(database)"
yamlKey="${10}"
# Example: "./scripts/test_plasmids.sh"
scriptToTest="${11}"
# Example: "/vol/spool/peter/meta-omics-toolkit/plasmid_yaml_database_tests"
workDirBasePath="${12}"
# Example: https,s3File
skipTests="${13:-noSkip}"
# Example: true
deleteScratchDir="${14:-no}"
# Example: --steps.export.emgb.kegg.database = "/path/to/db"
additionalYAMLParams="${15}"

mkdir -p ${basePath}

# Create all possible yaml configuration files
bash ./scripts/createDatabaseYML.sh "${extractedPath}" "${md5sum}" "${https}" \
	"${compressedPath}" "${s3File}" "${s3Directory}"  \
	"${s5cmdCommand}" \
	"${basePath}" "${yamlToTest}" "${yamlKey}" "${skipTests}"

# Run tests on all configuration files. The configurations are tested in parallel.
# Every run gets its own work directory, log file, output directory and database directory.
parallel="${MODULE_DB_TEST_PARALLEL:-4}"

resultDir="${workDirBasePath}/results"
rm -rf "${resultDir}"
mkdir -p "${resultDir}"

# Tests one configuration file. The output is written to a log file and the exit code to a separate file,
# so that all configurations are tested even if one of them fails.
function testConfiguration() {
	local yml="$1"
	local name
	name=$(basename "${yml}")
	local workDir="${workDirBasePath}/work_${name}"
	mkdir -p "${workDir}"
	local databaseDir="/vol/scratch/databases_$(cat /proc/sys/kernel/random/uuid)/"

	# Parallel runs must not resume the session of another run, and need their own nextflow log file and output directory.
	(
		export RESUME=""
		export NXF_LOG_FILE="${workDir}/nextflow.log"
		bash "${scriptToTest}" " --databases=${databaseDir} --output=${workDir}/output ${additionalYAMLParams} " "${yml}" "${workDir}" "slurm"
	) > "${resultDir}/${name}.log" 2>&1
	echo $? > "${resultDir}/${name}.exit"

	# Cleanup scratch directory
	if [[ "$deleteScratchDir" == "yes" ]]; then
		sinfo -o " %n %t"   | tail -n +2 | sed 's/^ //g' | grep -v " down" \
			| grep -v " drain" | grep -v "dummy" \
			| cut -d ' ' -f 1 | xargs -P 10 -I {}  ssh -tt ubuntu@{} " sudo rm -rf ${databaseDir} "
	fi
}
export -f testConfiguration
export resultDir workDirBasePath scriptToTest additionalYAMLParams deleteScratchDir

echo "Testing configuration files in ${basePath} (max. ${parallel} in parallel)"
find "${basePath}" -maxdepth 1 -type f | sort \
	| xargs -P "${parallel}" -I {} bash -c 'testConfiguration "$1"' _ {}

# Print the results. The logs of all configurations are printed as groups in GitHub Actions.
failed=0
for yml in $(find "${basePath}" -maxdepth 1 -type f | sort); do
	name=$(basename "${yml}")
	exitCode=$(cat "${resultDir}/${name}.exit" 2> /dev/null || echo "missing")
	if [[ "${GITHUB_ACTIONS}" == "true" ]]; then echo "::group::Configuration File ${name} (exit code: ${exitCode})"; fi
	echo "Configuration File ${yml}"
	cat "${resultDir}/${name}.log" 2> /dev/null
	if [[ "${GITHUB_ACTIONS}" == "true" ]]; then echo "::endgroup::"; fi
	if [[ "${exitCode}" != "0" ]]; then failed=$((failed + 1)); fi
done

echo "Summary:"
for yml in $(find "${basePath}" -maxdepth 1 -type f | sort); do
	name=$(basename "${yml}")
	exitCode=$(cat "${resultDir}/${name}.exit" 2> /dev/null || echo "missing")
	if [[ "${exitCode}" == "0" ]]; then echo "  PASSED ${name}"; else echo "  FAILED ${name} (exit code: ${exitCode})"; fi
done

if [[ "${failed}" -gt 0 ]]; then
	echo "${failed} configuration(s) failed."
	exit 1
fi
