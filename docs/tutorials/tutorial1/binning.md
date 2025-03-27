
We will now perform binning of our assembly.

Binning is a critical process in metagenomics that involves grouping DNA fragments, known as contigs,
into bins that likely originate from the same organism or genome.
This technique is essential for reconstructing genomes from mixed microbial communities,
where DNA from various organisms is intermingled,
and there is often no prior knowledge of the species present.

The binning process primarily utilizes two types of information: composition and coverage of the assembled contigs.

 * **Composition:** This refers to the analysis of k-mer frequency profiles within the contigs. K-mers are short DNA sequences of a fixed length (commonly tetranucleotides, or sequences of four nucleotides). These k-mer frequencies are generally conserved within a single genome, meaning that different regions of the same genome will share similar sequence composition.

 * **Coverage:** This reflects the abundance of each contig in the assembly. Organisms that are more abundant in the environment contribute more DNA to the metagenomic sample,
resulting in their sequences being represented more frequently in the data.

During binning, binning tools look for DNA fragments that, despite not being assembled together,
exhibit similar composition and occur at approximately equal abundances.
These similarities suggest that such contigs likely originate from the same genome,
enabling accurate grouping into bins.

By integrating both composition and coverage information, 
binning provides a powerful approach to distinguish between different organisms within a mixed sample. This method is particularly valuable as it allows researchers to reconstruct genomes without prior knowledge of the microbial community, making it an indispensable tool in metagenomic studies.


## MetaBAT

MetaBAT is an automated metagenome binning software
which integrates empirical probabilistic distances of genome abundance
and tetranucleotide frequency. See the [MetaBAT home page](https://bitbucket.org/berkeleylab/metabat) for more information.
  
We will now perform a binning of our assembly with MetaBAT using the Metagenomics-Toolkit.
The following snippet represents the Toolkit configuration for the binning module:

```YAML linenums="1" title="Binning Configuration File Snippet 1"
---8<--- "default/tutorials/tutorial1/fullpipeline_binning.yml:47:60"
```

!!! question "Task 1"
    We can run it directly with:
    ```BASH
    ---8<--- "scripts/tutorials/tutorial1/test_binning.sh:3:12"
    ```

For reference, your complete parameter file looks like this:
??? Parameter-file

    ```BASH
    ---8<--- "default/tutorials/tutorial1/fullpipeline_binning.yml"
    ```    

## Assembly Evaluation by Read Mapping

Before we take a closer look at the binning results, 
let us first check the mapping of the reads to the assembly (alignment) that was computed as part of the binning process and saved in the file with the `.bam` prefix.
The number of reads mapped to the contigs of the assembly is called the contig abundance and is also used as input to a binning tool to separate contigs.
From this alignment we get the contig **coverage** (read depth) which refers to the number of times each nucleotide position in a contig is covered by sequencing reads on average.
It quantifies how many reads align to each base pair in the assembled sequence.

!!! question "Task 2"
    First we have to index the sorted BAM file:
    ```BASH
    cd ~/mgcourse/
    ```
    ```BASH
    samtools index output/data/1/binning/*/contigMapping/data.bam
    ```

Usually a bam file has to be sorted before it can be indexed. The sorting has already been done for you by the Toolkit.
    
!!! question "Task 3"
    To look at the BAM file you can use the following command (You can quit less with `q`):
    ```BASH
    cd ~/mgcourse/
    ```
    ```BASH
    samtools view output/data/1/binning/*/contigMapping/data.bam | less
    ```

IGV, or Integrated Genome Viewer, is an open-source bioinformatics tool designed for the interactive visualization and exploration of genomic data.
It primarily serves as a platform to view alignments of sequencing reads against a reference genome,
aiding in the identification of genetic variants such as SNPs (Single Nucleotide Polymorphisms) and indels (insertions and deletions).
Beyond read alignment, IGV supports various data types, including gene annotations and expression levels, 
offering a comprehensive view that enhances understanding of the genomic context.
In our case, we are only interested in the alignment.
    
!!! question "Task 4"
    Now copy/link everything you need for igv in a separate folder:
    
    ```BASH
    cd ~/mgcourse/
    ```
    ```BASH
    mkdir igv_data
    cd igv_data
    cp ../output/data/1/assembly/*/megahit/data_contigs.fa.gz .
    gunzip data_contigs.fa.gz
    ln -s ../output/data/1/binning/*/contigMapping/data.bam
    ln -s ../output/data/1/binning/*/contigMapping/data.bam.bai
    ```
       
!!! question "Task 5"
    We will use the IGV genome browser to look at the mappings.
    ```BASH
    igv
    ```
    
    Now let's look at the mapped reads:
    
    1. Load the contig sequences into IGV. Use the menu `Genomes->Load Genome from File...`
    2. Load the BAM file into IGV. Use menu `File->Load from File...`
    
!!! question "Task 6"
    Look for errors in the mappings - are all those error sequencing errors?
    ??? Solution
        Some errors are due to merging multiple strains into one contig. This can be clearly seen 
        when there are a large number of errors at one position with different possibilities for that base.

## MetaBAT Results

Let's now inspect the bins created by MetaBAT:

!!! question "Task 7"
    How many bins did metabat generate? Locate the metabat results and the fasta files for each bin in the output folder.
    ??? Solution
        ```BASH
        cd ~/mgcourse/
        ```
        ```BASH
        ls -l output/data/1/binning/*/metabat/
        ```
        ```BASH
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
        MetaBAT generated 10 bins.

You might have noticed that there is also a `data_notBinned.fa` file that represents contigs that could not be binned. We will investigate these contigs in the
next section of this tutorial.
The next question you might want to ask is whether you can trust these bins and which organism they represent according to a taxonomy.
Before that, let's have a look at what other information the Metagenomics-Toolkit provides as part of the binning output.

!!! question "Task 8"
    There is a file containing some statistics on the generated Bins `data_bins_stats.tsv`. 
    Find out, which bin has the highest coverage and which one has the highest N50.
    ??? Solution    
        You can have a look on some bin statistics with:
        ```
        cd ~/mgcourse/
        ```
        ```BASH
        less output/data/1/binning/0.5.0/metabat/data_bins_stats.tsv
        ```
        ```BASH
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
        In this result, bin 8 has the highest coverage (15.5364) and highest N50 (13397). Also note, that some of the bins highly differ in their GC content.

---

➡️ [**Continue to: Assessing Bin Quality**](./bin_quality.md)

