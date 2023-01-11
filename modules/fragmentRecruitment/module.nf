include { wMashScreenFile; wMashScreenList; } from  './mashScreen'
include { wFragmentRecruitmentFile as wFragmentRecruitmentFileFrhit  } from './frhit'


workflow wFragmentRecruitmentFile {
   main:
     wMashScreenFile()
     if(params.steps.fragmentRecruitment.containsKey("frhit")){
       wFragmentRecruitmentFile(Channel.fromPath(params?.steps?.fragmentRecruitment?.frhit?.samples), Channel.fromPath(params?.steps?.fragmentRecruitment?.frhit?.genomes))
     }
}


workflow wFragmentRecruitmentList {
   take:
     pairedReads
     ontReads
     medianQuality
   main:
     wMashScreenList(pairedReads, ontReads, medianQuality)
   emit:
     binsStats = wMashScreenList.out.binsStats
     genomes = wMashScreenList.out.genomes
     foundGenomesSeperated = wMashScreenList.out.foundGenomesSeperated
     contigCoverage = wMashScreenList.out.contigCoverage
}
