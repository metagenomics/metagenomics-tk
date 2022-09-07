
process pDumpLogs {

    publishDir params.output, mode: "${params.publishDirMode}", saveAs: { filename -> "${outputDir}/" + filename }, \
        pattern: "{**.out,**.err,**.sh,**.log}"

    input:
    tuple env(ID), val(outputDir), val(maxLogLevel) ,file("command.sh"), file("command.out"), file("command.err"), file("command.log")

    output:
    path("logs/*"), optional: true, emit: logs

    shell:
    '''
    if [ !{maxLogLevel} -ge !{params.logLevel} ];
    then
       mkdir -p logs
       cp command.log  logs/${ID}.log
       cp command.err  logs/${ID}.err
       cp command.out  logs/${ID}.out
       cp command.sh  logs/${ID}.sh
    fi
    '''
}

