include { wMultiBinningShortReadList; } from './multiBinning'
include { wShortReadBinningList; _wPostProcessBins; } from './shortReadBinning'
include { wRefinementList } from './binRefinement.nf'
include { wSaveSettingsList } from '../config/module'


workflow wHydraBinFile {
    SAMPLE_IDX = 0
    SAMPLE_PAIRED_IDX = 1
    BINNING_LABEL_IDX = 1
    UNPAIRED_IDX = 2

    channel.from(file(params.steps.hydraBin.input.samples))
        | splitCsv(sep: '\t', header: true)
        | multiMap { row ->
            reads: [row.SAMPLE, row.READS, file('empty')]
            contigs: [row.SAMPLE, row.CONTIGS]
            binningLabels: [row.SAMPLE, row.MULTI_BINNING_GROUP]
            gfa: [row.SAMPLE, row.GFA]
            paths: [row.SAMPLE, row.PATHS]
        } | set { samples } 

    samples.binningLabels | groupTuple(by: BINNING_LABEL_IDX) 
         | flatMap { samples, group -> 
            numberOfSamples = samples.size()
            samples.collect { sample -> 
                [ sample, group, numberOfSamples > 1 ? true : false, numberOfSamples ] 
            }
          } | set { binningLabels }

    wSaveSettingsList(samples.reads | map { it -> it[SAMPLE_IDX] })

    _wProcessShortReadBinning(samples.reads, samples.contigs, binningLabels, samples.gfa, samples.paths)
}


workflow _wProcessShortReadBinning {
      take:
        qcReads
        contigs
        binningLabels
        gfa
        paths
      main:
      // Figure out whether the sample belongs to a multi binning group
      //
      DO_NOT_ESTIMATE_IDENTITY = "-1"
      SAMPLE_IDX=0
      GROUP_IDX=0
      GROUP_2_IDX=1
      GROUP_SIZE_IDX=3
      IS_MULTI_SAMPLE_IDX = 4
      READS_FILE_IDX = 1
      READS_UNPAIRED_FILE_IDX = 2
      qcReads | combine(binningLabels, by:SAMPLE_IDX ) | branch { sample ->
        singleSample: !sample[IS_MULTI_SAMPLE_IDX]
        multiSample: sample[IS_MULTI_SAMPLE_IDX]
      } | set {sampleTypeReads}

      // Make sure that the number of contigs per group matches the number of read samples per group because certain samples may fail.
      sampleTypeReads.multiSample 
        | map { sample -> [sample[SAMPLE_IDX], sample[READS_FILE_IDX], sample[READS_UNPAIRED_FILE_IDX]] } 
        | combine(contigs, by: SAMPLE_IDX) | set { qualityCheckedData }
      
      // Check the number of samples per group
      // If there is more than one sample then go on to multibinning
      // otherwise continue with single binning
      qualityCheckedData | join(binningLabels | map { sample -> [sample[SAMPLE_IDX], sample[GROUP_2_IDX], sample[GROUP_SIZE_IDX]]}, by: SAMPLE_IDX) 
        | map { sample, readsPair, readsSingle, contigs, group, groupSize -> tuple( groupKey(group, groupSize), [sample, readsPair, readsSingle, contigs, group, groupSize]) }
        | groupTuple(remainder: true)
        | branch { group, samples ->
            singleSample: samples.size() == 1
            multiSample: samples.size() > 1
            failedSample: samples.size() == 0 
        } | set { checkedSamples }

        checkedSamples.multiSample | map{ group, samples -> samples} | flatMap | set { multiSamples}

        checkedSamples.singleSample | map{ group, samples -> samples} | flatMap | set { singleSample}

        multiSamples | multiMap { sample, readsPair, readsSingle, contigs, group, groupSize  ->
            contigs: [sample, contigs]
            reads: [sample, readsPair, readsSingle]
            binningLabels: [sample, group, groupSize]
        } | set { multiSamplesInput } 

      wMultiBinningShortReadList(multiSamplesInput.contigs, multiSamplesInput.reads, multiSamplesInput.binningLabels)

      singleSample 
      | combine(gfa, by: SAMPLE_IDX)
      | combine(paths, by: SAMPLE_IDX)
      | multiMap { sample, readsPair, readsSingle, contigs, group, groupSize, gfa, paths ->
            contigs: [sample, contigs]
            gfa: [sample, gfa]
            paths: [sample, paths]
            reads: [sample, readsPair, readsSingle]
      } | set { singleSampleInput } 

      contigs 
      	| combine(binningLabels 
	      | filter({sample, group, isMultiSample, groupCount -> !isMultiSample }), by: SAMPLE_IDX)  
	      | map { sample, contigs, group, isMultiSample, groupCount -> [sample, contigs] }  
      	| set { singleSampleContigs }

      gfa
      	| combine(binningLabels 
	      | filter({sample, group, isMultiSample, groupCount -> !isMultiSample }), by: SAMPLE_IDX)  
	      | map { sample, gfa, group, isMultiSample, groupCount -> [sample, gfa] }  
      	| set { singleSampleGfa }

      paths
      	| combine(binningLabels 
	      | filter({sample, group, isMultiSample, groupCount -> !isMultiSample }), by: SAMPLE_IDX)  
	      | map { sample, paths, group, isMultiSample, groupCount -> [sample, paths] }  
      	| set { singleSamplePaths }

      sampleTypeReads.singleSample
            | map { sample -> [sample[SAMPLE_IDX], sample[READS_FILE_IDX], sample[READS_UNPAIRED_FILE_IDX]] } | mix(singleSampleInput.reads)
            | set { singleSampleReadInput }

      singleSampleGfa | mix(singleSampleInput.gfa) | set {singleSampleGfaInput}  

      singleSamplePaths | mix(singleSampleInput.paths) | set {singleSamplePathsInput}

      singleSampleContigs 
            | mix(singleSampleInput.contigs)
            | set { singleSampleContigsInput }

      binContigMapping = channel.empty()


      if(params.steps.containsKey("binRefinement") 
        && params.steps.binRefinement.mode.includeMultiSample){

        singleSampleReadInput 
            | mix(multiSamplesInput.reads)
            | set { singleSampleReadInput }

        contigs 
		| set { singleSampleContigsInput }

        multiSamples
          | combine(gfa, by: SAMPLE_IDX)
          | combine(paths, by: SAMPLE_IDX)
          | multiMap { sample, readsPair, readsSingle, contigs, group, groupSize, gfa, paths ->
                contigs: [sample, contigs]
                gfa: [sample, gfa]
                paths: [sample, paths]
                reads: [sample, readsPair, readsSingle]
        } | set { singleSampleShorReadInput } 

       singleSampleGfaInput | mix(singleSampleShorReadInput.gfa) 
           | set { singleSampleGfaInput }

       singleSamplePathsInput | mix(singleSampleShorReadInput.paths)
           | set { singleSamplePathsInput }

       binContigMapping | mix(wMultiBinningShortReadList.out.binContigMapping)
           | set { binContigMapping } 
      } 

      wShortReadBinningList(singleSampleContigsInput,  
        singleSampleReadInput)

      wShortReadBinningList.out.binsStatsInput | set { binsStatsInputShort }

      wShortReadBinningList.out.mapping | set { mappingShort }

      wShortReadBinningList.out.bins 
	| mix(wMultiBinningShortReadList.out.bins) | set { bins }

      wShortReadBinningList.out.notBinnedContigs  | set {notBinnedContigs}

      binsStats = channel.empty()

      if (params.steps.containsKey("binRefinement")) {
        wRefinementList(singleSampleContigsInput, 
        binContigMapping | mix(wShortReadBinningList.out.binContigMapping),
        bins,
        binningLabels,
        notBinnedContigs,
        singleSampleGfaInput,
        singleSamplePathsInput,
        qcReads) 
        
        wRefinementList.out.bins | set { bins }
        wRefinementList.out.notBinned | set { notBinnedContigs }

        wRefinementList.out.binContigMapping
            | map { sample, method, binContigMapping -> [sample, binContigMapping]}
            | join(mappingShort, by: SAMPLE_IDX)
            | combine(channel.from("refinement/final"))
            | join(bins | map { sample, method, bins -> [sample, bins]}, by: SAMPLE_IDX)
            | combine(channel.value(DO_NOT_ESTIMATE_IDENTITY))
            | set { binsStatsInputShort }
      }  

      if(!params.steps.containsKey("binRefinement") || 
        !params.steps.binRefinement.mode.includeMultiSample){
                
              binsStats | mix(wMultiBinningShortReadList.out.binsStats) 
                | set {binsStats}

              bins | mix(wMultiBinningShortReadList.out.bins)
                | set { bins }

              notBinnedContigs | mix(wMultiBinningShortReadList.out.notBinnedContigs)
                | set { notBinnedContigs }
      }

      _wPostProcessBins(mappingShort, binsStatsInputShort, 
      bins | map { sample, method, bins -> [sample, bins]})
      _wPostProcessBins.out.binMap | mix(binsStats) | set { binsStats }

      bins | map { sample, method, bins -> [sample, bins]} 
        | set { bins }

      notBinnedContigs | map { sample, method, notBinned -> [sample, notBinned]} 
        | set { notBinnedContigs }

      mappingShort | mix(wMultiBinningShortReadList.out.mapping)
        | set { mapping }

      wShortReadBinningList.out.unmappedReads
        | mix(wMultiBinningShortReadList.out.unmappedReads)
        | set { unmappedReads }

      wShortReadBinningList.out.contigCoverage
        | mix(wMultiBinningShortReadList.out.contigCoverage)
        | set { contigCoverage }

    emit:
      bins = bins 
      binsStats = binsStats
      notBinnedContigs = notBinnedContigs
      mapping = mapping
      unmappedReads = unmappedReads
      contigCoverage = contigCoverage
}
