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

              if(config.download.containsKey("s5cmd") && config.download.s5cmd.containsKey("keyfile")){
                  volumeMountStr += " --volume " + config.download.s5cmd.keyfile + ":/.aws/credentials" + " --volume " + config.download.s5cmd.keyfile + ":/root/.aws/credentials"
              }

              return volumeMountStr;
          }

      } else {
          return "";
      }
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
