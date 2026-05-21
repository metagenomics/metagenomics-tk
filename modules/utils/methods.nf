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


