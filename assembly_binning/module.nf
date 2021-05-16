nextflow.enable.dsl=2

params.output = "out"
params.input = "sample_all.tsv"
params.quast = false
params.spades = false
params.metaspades = false
params.megahit = false
params.checkm = false
params.metabat = false
params.maxbin = false
params.getreads = false
params.sra = false
params.interleaved = false
params.deinterleaved = false
params.gtdb_database = ""
params.checkm_database = ""
params.gtdb = false
params.skip = false
params.postprocessing.buffer = 30
params.ending = ".fa"


include { run_postprocess } from './postprocessing/module.nf'
include { run_binning } from './binning/module.nf'
include { runMegahitInterleaved; runBBMapDeinterleave; run_assembly;runMegahitSplit; } from './assembly/module.nf'



workflow assembly_binning_input {
     take:
       reads
     main:
       if(params.mode.interleaved){
          if(params.mode.skip){
            reads | splitCsv(sep: '\t', header: true) \
             | map { it -> [ it.SAMPLE, it.READS ]} \
             | runMegahitInterleaved     

            runMegahitInterleaved.out.reads_processed | set { processed_reads}
            run_binning(runMegahitInterleaved.out.contigs, runMegahitInterleaved.out.reads_processed)
            run_postprocess(run_binning.out.bins, run_binning.out.mapping)

            runMegahitInterleaved.out.fastp_summary | collectFile(newLine: false, keepHeader: true, storeDir: params.output + "/summary/" ){ item ->
              [ "fastp_summary.tsv", item[1].text  ]
            }
          } else {
            reads | splitCsv(sep: '\t', header: true) \
             | map { it -> [ it.SAMPLE, it.READS ]} \
             | runBBMapDeinterleave | run_assembly

            run_binning(run_assembly.out.contigs, run_assembly.out.reads_processed)
            run_postprocess(run_binning.out.bins, run_binning.out.mapping)

            run_assembly.out.fastp_summary | collectFile(newLine: false, keepHeader: true, storeDir: params.output + "/summary/"  ){ item ->
              [ "fastp_summary.tsv", item[1].text  ]
            }

          }
       }

       if(!params.mode.interleaved){
          if(params.mode.skip){

            reads | splitCsv(sep: '\t', header: true) \
             | map { it -> [ it.SAMPLE, it.READS1, it.READS2 ]} \
             | runMegahitSplit 
             
            runMegahitSplit.out.reads_processed | set { processed_reads}

            run_binning(runMegahitSplit.out.contigs, runMegahitSplit.out.reads_processed)
            run_postprocess(run_binning.out.bins, run_binning.out.mapping)
            runMegahitSplit.out.fastp_summary | collectFile(newLine: false, keepHeader: true, storeDir: params.output + "/summary/" ){ item ->
              [ "fastp_summary.tsv", item[1].text ]
            }
          } else {
            reads | splitCsv(sep: '\t', header: true) \
             | map { it -> [ it.SAMPLE, it.READS1, it.READS2 ]} \
             | run_assembly

            run_binning(run_assembly.out.contigs, run_assembly.out.reads_processed)
            run_postprocess(run_binning.out.bins, run_binning.out.mapping)

            run_assembly.out.fastp_summary | collectFile(newLine: false, keepHeader: true, storeDir: params.output + "/summary/"  ){ item ->
              [ "fastp_summary.tsv", item[1].text  ]
            }
          }
       }
    emit:
      bins = run_postprocess.out.bins_info
      processed_reads = processed_reads
}


workflow assembly_binning_input_sra  {
     if(params.sra){
        Channel.fromSRA(['SRR6820513'], apiKey: '')  | set{ input_reads }
     }
     //assembly_binning(input_reads)
}
