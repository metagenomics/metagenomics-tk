class Utils {


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
}
