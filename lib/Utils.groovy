class Utils {

  /**
  *
  * Get Docker mount point string for database folder if the file must be downloaded first,
  * otherwise mount the file directly if it is already available on the filesystem.
  *
  **/
  static String getDockerMount(config, params, useParentDirectory=false) {
      def getPathWithoutFile = { filePath ->
          int dotIndex = filePath.lastIndexOf('/');
          return (dotIndex == -1) ? filePath : filePath.substring(0, dotIndex);
      }

      def selectInput = { f -> useParentDirectory ? getPathWithoutFile(f): f }
      if(config!=null){
          if(config.containsKey("extractedDBPath")){
              return " --volume " + selectInput(config.extractedDBPath) + ":" + selectInput(config.extractedDBPath) ;
          } else if (config.containsKey("download")) {
              def volumeMountStr = ""
              if(config.download.source.startsWith("/")){
                  volumeMountStr += " --volume " +  selectInput(config.download.source) + ":" + selectInput(config.download.source)
              }
              volumeMountStr += " --volume " + params.polished.databases + ":" + params.polished.databases ;

              return volumeMountStr;
          }

      } else {
          return "";
      }
  }

  static String getCreateDatabaseDirCommand(db){
    return "if [ ! -z " + db + " ]; then mkdir " + db + " -p; fi"
  }


  static String getModulePath(module){
    return module.name + '/' + module.version.major + "." +
          module.version.minor + "." +
          module.version.patch
  }


  static String getBeforeScript(script, image){
    if(script.isEmpty()){
      return "echo 'No BeforeScript'";
    } else {
      return "bash " + script + " " + image ; 
    }
  }


  static Collection asList(element){
     if(element instanceof Collection){
         return element;
     } else {
         return [element];
     }
  }

  static String getPathWithoutFile(String filePath) {
    int dotIndex = filePath.lastIndexOf('/');
    return (dotIndex == -1) ? filePath : filePath.substring(0, dotIndex);
  }


  /*
   * This function sets a time limit based on a user provided mode and the resource defaults
   * that are usually consumed by the process.
   * 
   * Input:
   * process: Defines the user provided timeLimit mode.
   * processDefaults: Defines default resource values (cpus, RAM, etc.).
   * flavor: Defines the user provided flavor for the process. 
  */
  static String setTimeLimit(process, processDefaults, flavor){
    def timeUnit = "h"
    switch(process.timeLimit){
      case "DISABLED":
        return "";
      case "AUTO":
        def defaultCPUs = processDefaults.flavor.cpus
        if(defaultCPUs > flavor.cpus){
          def timeLimit = defaultCPUs/flavor.cpus * processDefaults.time;
          return timeLimit + timeUnit;
        } else {
          return processDefaults.time + timeUnit;
        }
      default:
        if(process.timeLimit instanceof Number){
          return process.timeLimit + timeUnit;
        } else {
          System.out.println("WARNING: unknown time limit parameter specified!");
        }
    }
  }

  static Object[] flattenTuple(tupl){
  	def chunkList = [];
  	def SAMPLE_IDX = 0;
  	def PATHS_IDX = 1;
  	tupl[PATHS_IDX].each {
     		chunkList.add([tupl[SAMPLE_IDX], it]);
  	}
  	return chunkList;
  }

  /*
  * This method takes number of entries in a input file (e.g. fata entries in multi fasta file),
  * the maximum number of allowed entries per chunk and the actual input (e.g. file).
  * It creates a list of indices of chunks of the input file based on the input parameters.
  */
  static List splitFilesIndex(seqCount, chunkSize, sample){
    int chunk=seqCount.intdiv(chunkSize)
    if(seqCount.mod(chunkSize) != 0){
      chunk = chunk + 1
    }
    def chunks = []
    for(def n : 1..chunk){
      int start = (n-1) * chunkSize + 1
     
      int end = n * chunkSize
    
      if(end > seqCount){
          end=seqCount
      }
      chunks.add(sample + [start, end])
    }

    return chunks
  }


  static getMappingIdentityParam(medianQuality) {
    if(medianQuality > 17){
      return 97
    }
    if(medianQuality > 13){
      return 95
    } else {
      return 90
    }
  }
}
