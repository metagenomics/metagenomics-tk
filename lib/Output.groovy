class Output {
  static String getOutput(sample, runid, tool, module, filename){
    return sample + '/' + runid + '/' + module.name + '/' + 
           module.version.major + "." +
           module.version.minor + "." +
           module.version.patch + 
           '/' + tool + '/' + filename
  }

}
