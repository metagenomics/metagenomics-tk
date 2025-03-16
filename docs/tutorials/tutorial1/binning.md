We will now run a binning of our assembly.

##Metabat

MetaBAT, An Efficient Tool for Accurately Reconstructing Single
Genomes from Complex Microbial Communities.

Grouping large genomic fragments assembled from shotgun metagenomic
sequences to deconvolute complex microbial communities, or metagenome
binning, enables the study of individual organisms and their
interactions. MetaBAT is an automated metagenome binning software
which integrates empirical probabilistic distances of genome abundance
and tetranucleotide frequency. See the [MetaBAT home page](https://bitbucket.org/berkeleylab/metabat>) for more info.
  
We will now run a binnning using metabat, these lines need to be added below the assembly part of the parameter file:

```BASH
---8<--- "default/tutorials/tutorial1/fullpipeline_binning.yml:41:54"
```

The complete parameter file is appended at the end of this page.

!!! question "Task 1"
    We can run it directly with:
    ```BASH
    ---8<--- "scripts/tutorials/tutorial1/test_binning.sh:3:12"
    ```

## Assembly evaluation by read mapping

Before we have a look on the binning results, we will have a look on the assembly results again, we now have a mapping of the reads back to the assembly which has been computed for the binning process. We will load that data after some preprocessing into the IGV genome browser.

!!! question "Task 2"
    We first have to sort the BAM file by starting position of the alignments. This can be done using samtools.
    ```BASH
      cd ~/mgcourse/
      samtools sort -o output/data/1/binning/0.5.0/contigMapping/data_sorted.bam -@ 28 output/data/1/binning/0.5.0/contigMapping/data.bam 
    ```

!!! question "Task 3"
    Now we have to index the sorted BAM file:
    ```BASH
      cd ~/mgcourse/
      samtools index output/data/1/binning/0.5.0/contigMapping/data_sorted.bam
    ```
    
!!! question "Task 4"
    To look at the BAM file use:
    ```BASH
      cd ~/mgcourse/
      samtools view output/data/1/binning/0.5.0/contigMapping/data_sorted.bam | less
    ```
    
!!! question "Task 5"
    Now copy/link everything you need for igv in a seperate folder:
    
    ```BASH
      cd ~/mgcourse/
      mkdir igv_data
      cd igv_data
      cp ../output/data/1/assembly/1.2.1/megahit/data_contigs.fa.gz .
      gunzip data_contigs.fa.gz
      ln -s ../output/data/1/binning/0.5.0/contigMapping/data_sorted.bam
      ln -s ../output/data/1/binning/0.5.0/contigMapping/data_sorted.bam.bai
    ```
       
!!! question "Task 6"
    We will use the IGV genome browser to look at the mappings.
    ```BASH
      igv
    ```
    
    Now let's look at the mapped reads:
    
    1. Load the contig sequences into IGV. Use the menu `Genomes->Load Genome from File...`
    2. Load the BAM file into IGV. Use menu `File->Load from File...`
    
!!! question "Task 7"
    Look for errors in the mappings - are all those error sequencing errors?
    ??? Solution
        Some errors are due to merging multiple strains into one contig, this can be clearly seen, when there are a large number of errors at one position with different possibilities for that base.

## Metabat results

MetaBAT generated genome BINs for our assembly. 

!!! question "Task 8"
    How many BINs did metabat generate? Locate the metabat results and the fasta files for each bin in the output folder.
    ??? Solution
        ```BASH
        cd ~/mgcourse/
        ls -l output/data/1/binning/0.5.0/metabat/
        data_bin.1.fa
        data_bin.10.fa
        data_bin.2.fa
        data_bin.3.fa
        data_bin.4.fa
        data_bin.5.fa
        data_bin.6.fa
        data_bin.7.fa
        data_bin.8.fa
        data_bin.9.fa
        data_bin_contig_mapping.tsv
        data_bins_stats.tsv
        data_contigs_depth.tsv
        data_notBinned.fa
        ```
        So MetaBAT did generate 10 contigs.

!!! question "Task 8"
    There is also a file containing some statistics on the generated BINs `data_bins_stats.tsv`. Find out, which BIN has the highest coverage and which hast the highest N50.
    ??? Solution    
        You can have a look on some bin statistics with:
        ```
        cd ~/mgcourse/
        less output/data/1/binning/0.5.0/metabat/data_bins_stats.tsv
        SAMPLE  BIN_ID  format  type    num_seqs        sum_len min_len avg_len max_len Q1      Q2      Q3      sum_gap N50     Q20(%)  Q30(%)  GC(%)   COVERAGE
        data    data_bin.1.fa   FASTA   DNA     66      682946  2837    10347.7 31045   6083.0  8277.5  13778.0 0       13397   0.00    0.00    25.03   15.5969
        data    data_bin.10.fa  FASTA   DNA     280     944933  2505    3374.8  6899    2786.5  3165.0  3690.5  0       3346    0.00    0.00    57.96   6.37826
        data    data_bin.2.fa   FASTA   DNA     150     520975  2500    3473.2  9084    2739.0  3127.0  3909.0  0       3328    0.00    0.00    43.93   7.13884
        data    data_bin.3.fa   FASTA   DNA     182     654974  2504    3598.8  8372    2821.0  3173.5  4021.0  0       3502    0.00    0.00    49.41   7.81954
        data    data_bin.4.fa   FASTA   DNA     139     1149027 2523    8266.4  35033   4367.0  6419.0  10607.0 0       10241   0.00    0.00    40.52   14.9085
        data    data_bin.5.fa   FASTA   DNA     454     1994574 2501    4393.3  13793   2998.0  3736.5  5327.0  0       4624    0.00    0.00    27.97   9.95115
        data    data_bin.6.fa   FASTA   DNA     88      280467  2509    3187.1  5229    2746.5  2998.5  3429.5  0       3093    0.00    0.00    39.02   5.18867
        data    data_bin.7.fa   FASTA   DNA     153     517295  2504    3381.0  7616    2769.0  3134.0  3688.0  0       3294    0.00    0.00    33.15   6.26526
        data    data_bin.8.fa   FASTA   DNA     180     958807  2514    5326.7  17160   3436.0  4467.0  6472.5  0       5965    0.00    0.00    33.35   10.8953
        data    data_bin.9.fa   FASTA   DNA     224     717735  2500    3204.2  6592    2692.0  2958.5  3442.0  0       3109    0.00    0.00    50.72   5.11192
        ```
        In this result, bin 1 has the highest coverage (15.5969) and bin 4 hast the highes N50 (10241). Also note, that some of the bins highly differ in their GC content.


        
For reference, your complete parameter file should look like this:
??? Parameter-file

    ```BASH
    ---8<--- "default/tutorials/tutorial1/fullpipeline_binning.yml"
    ```       

