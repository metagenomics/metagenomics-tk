class Utils {

  /**
  *
  * Get Docker mount point string for database folder if the file must be downloaded first,
  * otherwise mount the file directly if it is already available on the filesystem.
  *
  **/
  static String getDockerMount(config, params) {

    if(config!=null){
        if(config.containsKey("extractedDBPath")){
                return " --volume " + config.extractedDBPath + ":" + config.extractedDBPath ;
        } else if (config.containsKey("download")) {
                def volumeMountStr = ""
                if(config.download.source.startsWith("/")){
                        volumeMountStr += " --volume " + config.download.source + ":" + config.download.source  
                }
                volumeMountStr += " --volume " + params.polished.databases + ":" + params.polished.databases ;

        	if(config.download.containsKey("s5cmd") && config.download.s5cmd.containsKey("keyfile")){
                	volumeMountStr += " --volume " + config.download.s5cmd.keyfile + ":/.aws/credentials   "
        	}

        	return volumeMountStr;
        }

    } else {
    	return "";
    }
  }

  static String getDockerMountMMseqs(config, params) {

    if(config!=null){
        if(config.containsKey("extractedDBPath")){ 
                return " --volume " + getPathWithoutFile(config.extractedDBPath) + ":" + getPathWithoutFile(config.extractedDBPath) ;
        } else if (config.containsKey("download")) {
                def volumeMountStr = ""
                if(config.download.source.startsWith("/")){
                        volumeMountStr += " --volume " + getPathWithoutFile(config.download.source) + ":" + getPathWithoutFile(config.download.source)
                }
                volumeMountStr += " --volume " + params.polished.databases + ":" + params.polished.databases ;

                if(config.download.containsKey("s5cmd") && config.download.s5cmd.containsKey("keyfile")){
                        volumeMountStr += " --volume " + config.download.s5cmd.keyfile + ":/.aws/credentials   "
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

}
