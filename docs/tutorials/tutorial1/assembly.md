Once your raw reads have been quality-controlled, we will now perform an assembly using the Metagenomics-Toolkit.
The Toolkit supports several assemblers, one of which is the MEGAHIT assembler, which we will run in this section.

## MEGAHIT 

MEGAHIT is a single-node assembler for large and complex metagenomics
NGS reads, such as soil. It makes use of a succinct de Bruijn graph
(SdBG) to achieve a low memory assembly. See the
[MEGAHIT home page](https://github.com/voutcn/megahit/) for more
info.


## Metagenomics-Toolkit  

The following lines represent the part of the configuration that tells the Toolkit to run the assembly:

```YAML linenums="1" title="Assembly Configuration File Snippet 1"
---8<--- "default/tutorials/tutorial1/fullpipeline_assembly.yml:36:46"
```
!!! question "Task 1"

    Which tool is used for the assembly? What additional parameters are used?

    ??? Solution 
        MEGAHIT with the additional parameters: minimum contig length of `1000 bp` and the preset `meta-sensitive`.
        `meta-sensitive` sets the following parameters: `--min-count 1 --k-list 21,29,39,49,...,129,141`, which causes 
        MEGAHIT to use a longer list of k-mers.
    

!!! question "Task 2"

    Copy the following command to run the quality control and assembly module of the Toolkit on your machine
    
    ```BASH
    ---8<--- "scripts/tutorials/tutorial1/test_assembly.sh:3:13"
    ```

Note that you have used the `resume` flag this time, which causes Nextflow to reuse the results of the QC analysis. 
That's why you'll see several processes labeled **Cache process** on the Toolkit execution screen.

For reference, your complete parameter file looks like this:
??? Parameter-file

    ```BASH
    ---8<--- "default/tutorials/tutorial1/fullpipeline_assembly.yml"
    ```    

## Assembly Results

We will now have a first look at some assembly statistics. First of all, locate your assembly results somewhere in your `output` directory.

!!! question "Task 3"
    Where are the assembly results of MEGAHIT stored? And what files are generated?
    
    ??? Solution 
        The assembly results are stored in `output/data/1/assembly/1.2.1/megahit/`
        ```
        cd ~/mgcourse/
        ```
        ```BASH
        ls -l output/data/1/assembly/1.2.1/megahit/
        ```
        ```BASH
            data_contigs.fa.gz  # the assembled sequences as gzipped fasta
            data_contigs.fastg  # the assembly graph for inspection for example with Bandage
            data_contigs_stats.tsv  # some assembly statistics
        ```
!!! question "Task 4"
    Have a look at the `data_contigs_stats.tsv` file - how large is the assembly and what is the N50?
    
    ??? Solution
        ```BASH
        cd ~/mgcourse/
        ```
        ```BASH
        cat output/data/1/assembly/1.2.1/megahit/data_contigs_stats.tsv | column -s$'\t' -t
        ```
        ```BASH
        SAMPLE  file                format  type  num_seqs  sum_len   min_len  avg_len  max_len  Q1      Q2      Q3      sum_gap  N50   Q20(%)  Q30(%)  GC(%)
        data    data_contigs.fa.gz  FASTA   DNA   9172      20146206  1000     2196.5   35033    1274.0  1670.0  2450.0  0        2346  0.00    0.00    42.26
        ```
        In this assembly run, the assembly length is 20,146,206 bp and the N50 is 2,346 bp. 
 
        
---

➡️ [**Continue to: Assembly Evaluation**](./assembly_evaluation.md) 
    
