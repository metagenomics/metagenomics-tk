/*
*
* Method takes a list of lists of the form [[SAMPLE, BIN 1 path]] 
* and produces a map of the form [BIN_ID:bin.name, SAMPLE:sample, PATH:bin]
*/
def createMap(binning) {
    def chunkList = []
    binning.each {
        def sample = it[0]
        def bin = file(it[1])
        def binMap = [BIN_ID: bin.name, SAMPLE: sample, PATH: bin]
        chunkList.add(binMap)
    }
    return chunkList
}

/*
*
* Method takes two channels with map entries and two keys as input.
* Channels are joined by the keys provided.
* Resulting channel is returned as output.
*
*/
def mapJoin(channel_a, channel_b, key_a, key_b) {
    channel_a
        | map { it -> [it[key_a], it] }
        | cross(channel_b | map { it -> [it[key_b], it] })
        | map { it[0][1] + it[1][1] }
}

/*
*
* This function decides if the MMSeqs taxonomy process should be executed on MAGs.
*
*/
def runMMSeqsTaxonomy(isMMSeqsTaxonomySettingSet, isBinned, isPlasmid, runOnBinned) {
    // True if the MMSeqsTaxonomy setting is set and plasmids are not used as input and 
    // either the sample is not binned 
    // or should explicitly run on binned samples 
    return isMMSeqsTaxonomySettingSet && !isPlasmid && (!isBinned || runOnBinned)
}

/*
*
* The annotation step searches the proteins either with MMseqs2 or with DIAMOND, but never with both.
* Results of both tools would be mixed up in the steps that use them (e.g. KEGG and the export).
*
*/
def checkAnnotationSearchTool(annotationParams) {
    if (annotationParams.containsKey("mmseqs2") && annotationParams.containsKey("diamond")) {
        error "Only one of 'mmseqs2' and 'diamond' can be configured in the annotation step."
    }
    if (annotationParams.containsKey("keggFromMMseqs2")) {
        error "The annotation key 'keggFromMMseqs2' was renamed to 'keggFromBlast'."
    }
}
