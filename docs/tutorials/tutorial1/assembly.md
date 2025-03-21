We are going to run an assembly with the metagenomics toolkit using MEGAHIT as an assembler, the toolkit is also capable of running metaSPades. 

## MEGAHIT 

MEGAHIT is a single node assembler for large and complex metagenomics
NGS reads, such as soil. It makes use of succinct de Bruijn graph
(SdBG) to achieve low memory assembly. MEGAHIT can optionally utilize
a CUDA-enabled GPU to accelerate its SdBG contstruction. See the
[MEGAHIT home page](https://github.com/voutcn/megahit/) for more
info.


## Metagenomics-Toolkit  


The following lines need to be added to your parameter file in order to run the assembly (below the qc-part):
```BASH
---8<--- "default/tutorials/tutorial1/fullpipeline_assembly.yml:30:40"
```
!!! question "Task 1"

    Which tool is used for the Assembly? What additional parameters are used?

    ??? Solution 
        MEGAHIT with the additional parameters: minimum contig length of 1000bp and --meta-sensitive: which sets the following parameters: '--min-count 1 --k-list 21,29,39,49,...,129,141'
    

!!! question "Task 2"

    Copy the following command to run the quality control and assembly module of the Toolkit on one machine
    


    ```BASH
    ---8<--- "scripts/tutorials/tutorial1/test_assembly.sh:3:12"
    ```

## Assembly results
We will now have a first look on some assembly statistics. First of all, locate your assembly results somewhere in your `output` directory.

!!! question "Task 3"
    Where are the assembly results of MEGAHIT stored? And what files are generated?
    
    ??? Solution 
        The assembly results are stored in `output/data/1/assembly/1.2.1/megahit/`
        ```
        cd ~/mgcourse/
        ls -l output/data/1/assembly/1.2.1/megahit/
            data_contigs.fa.gz  # the assembled sequences as gzipped fasta
            data_contigs.fastg  # the assembly graph for inspection for example with Bandage
            data_contigs_stats.tsv  # some assembly statistics
        ```
!!! question "Task 4"
    Have a look at the `data_contigs_stats.tsv` file - how large is the assembly and what is the N50?
    
    ??? Solution
        ```
        cd ~/mgcourse/
        less output/data/1/assembly/1.2.1/megahit/data_contigs_stats.tsv
        
        SAMPLE  file    format  type    num_seqs        sum_len min_len avg_len max_len Q1      Q2      Q3      sum_gap N50     Q20(%)  Q30(%)  GC(%)
        data    data_contigs.fa.gz      FASTA   DNA     9172    20146206        1000    2196.5  35033   1274.0  1670.0  2450.0  0       2346    0.00    0.00    42.26
        ```
        In this assembly run, the assembly length is 20,146,206 bp and the N50 is 2,345 bp.
 
 
        
For reference, your complete parameter file should look like this:
??? Parameter-file

    ```BASH
    ---8<--- "default/tutorials/tutorial1/fullpipeline_assembly.yml"
    ```       
    
