/*
* This method collects files of the modules specifified by the "modules" parameter.  
* `dir` is the the path to the sra id
* `sra` is the SRA ID
*/
def collectModuleFiles(dir, sra, modules){
   def fileList = [];
   def moduleList = [];
   def moduleName = "";

   params.modules.eachWithIndex { v, k -> moduleList.add(v.getKey() + "/" + v.getValue().version.major + ".") }

   // iterate  over all specified modules
   for(module in modules){
     def moduleDir = file(dir + "/" + module.name + "/")
     moduleName = module.name
     // Check if the module exists
     if(moduleDir.exists()){
       // collect all files
       moduleDir.eachFileRecurse { item ->
           // make sure that only the module outputs of the specified version are collected.
           def found = moduleList.any {  item ==~ '.*' +  it + '.*'  }
           if(found){
              fileList.add([sra, item]);
           }
       }
     }
   }
   
   return fileList;
}

process pDumpLogs {

    tag "ID: $ID, Output: $outputDir, LogLevel: $maxLogLevel"

    publishDir params.output, mode: "${params.publishDirMode}", saveAs: { filename -> "${outputDir}/" + filename }, \
        pattern: "{**.out,**.err,**.sh,**.log}"

    input:
    tuple val(ID), val(outputDir), val(maxLogLevel) ,file("command.sh"), file("command.out"), file("command.err"), file("command.log")

    output:
    path("logs/*"), optional: true, emit: logs

    shell:
    '''
    if [ !{maxLogLevel} -ge !{params.logLevel} ];
    then
       mkdir -p logs
       cp command.log  logs/!{ID}.log
       cp command.err  logs/!{ID}.err
       cp command.out  logs/!{ID}.out
       cp command.sh  logs/!{ID}.sh
    fi
    '''
}


process pPublish {

    label 'tiny'

    container "${params.ubuntu_image}"

    publishDir "${outDir}", mode: "${params.publishDirMode}"

    errorStrategy 'retry'

    when:
    params?.input.containsKey("SRA")

    input:
    val(outDir)
    path(in)

    output:
    path("$in", includeInputs: true)

    shell:
    '''
    '''
}

