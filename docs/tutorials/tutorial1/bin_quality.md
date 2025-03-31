To get an impression of the quality of our bins, we compute the completeness and contamination values for our bins. 

## Computing Completeness and Contamination using CheckM2

CheckM 1 and 2 provide a set of tools to assess the quality of genomes recovered from isolates, single cells, or metagenomes.
While the first version relies on collocated sets of genes that are ubiquitous and single-copy within a phylogenetic lineage,
the second version uses a machine-learning-based approach that was trained on all available genomic information such as multi-copy genes,
biological pathways and modules. Both versions are supported by the Metagenomics-Toolkit and we will use the second version for this part of the tutorial.


The following lines represent the part of the configuration that tells the Toolkit to run the MagAttributes module that also includes
the CheckM2 tool:

```YAML linenums="1" title="MagAttributes Configuration File Snippet 1"
---8<--- "default/tutorials/tutorial1/fullpipeline_bin_quality.yml:61:67"
```
    
For reference, your complete parameter file looks like this:
??? Parameter-file

    ```BASH
    ---8<--- "default/tutorials/tutorial1/fullpipeline_bin_quality.yml"
    ```    

!!! Question "Task 2"
    Run CheckM2 with the Metagenomics-Toolkit using the following command:
    ```BASH
    ---8<--- "scripts/tutorials/tutorial1/test_bin_quality.sh:3:14"
    ```


!!! Question "Task 3"
    Locate the CheckM results inside the `output` directory and find out the completeness and contamination values for all of our bins.
    ??? Solution
        ```BASH
        cd ~/mgcourse
        ```
        ```BASH
        cut -f2,3,4 output/data/1/magAttributes/*/checkm2/data_checkm2_generated.tsv
        ```
        ```BASH        
        BIN_ID	COMPLETENESS	CONTAMINATION
        data_bin.1.fa	16.88	0.0
        data_bin.10.fa	24.04	0.04
        data_bin.2.fa	12.81	0.12
        data_bin.3.fa	20.57	0.01
        data_bin.4.fa	52.73	0.3
        data_bin.5.fa	11.09	0.0
        data_bin.6.fa	81.91	0.46
        data_bin.7.fa	20.68	0.75
        data_bin.8.fa	85.4	0.28
        data_bin.9.fa	39.66	0.62
        ```        

!!! Question "Task 4"
    Compare the CheckM2 completeness results with the genome fraction results from QUAST. Do they match your expectations?
    ??? Solution
        No, completeness is much lower than expected in comparison to the genome fraction we observed in the QUAST results (see Assembly evaluation part).

We will now run QUAST again, but this time we specify the bins and the unbinned contigs instead of the assemblies as input.
By doing this, QUAST can tell us which reference genome the bin belongs to.

!!! Question "Task 5"
    Run metaquast with the individual bins as input:
    
    ```BASH
    cd ~/mgcourse/
    ```

    ```BASH
    metaquast.py --threads 28 --gene-finding \
    -R genomes/Aquifex_aeolicus_VF5.fna,\
    genomes/Bdellovibrio_bacteriovorus_HD100.fna,\
    genomes/Chlamydia_psittaci_MN.fna,\
    genomes/Chlamydophila_pneumoniae_CWL029.fna,\
    genomes/Chlamydophila_pneumoniae_J138.fna,\
    genomes/Chlamydophila_pneumoniae_LPCoLN.fna,\
    genomes/Chlamydophila_pneumoniae_TW_183.fna,\
    genomes/Chlamydophila_psittaci_C19_98.fna,\
    genomes/Finegoldia_magna_ATCC_29328.fna,\
    genomes/Fusobacterium_nucleatum_ATCC_25586.fna,\
    genomes/Helicobacter_pylori_26695.fna,\
    genomes/Lawsonia_intracellularis_PHE_MN1_00.fna,\
    genomes/Mycobacterium_leprae_TN.fna,\
    genomes/Porphyromonas_gingivalis_W83.fna,\
    genomes/Wigglesworthia_glossinidia.fna \
    -o quast_bins \
    -l BIN1,BIN2,BIN3,BIN4,BIN5,BIN6,BIN7,BIN8,BIN9,BIN10,UNBINNED \
    output/data/1/binning/*/metabat/data_bin.1.fa \
    output/data/1/binning/*/metabat/data_bin.2.fa \
    output/data/1/binning/*/metabat/data_bin.3.fa \
    output/data/1/binning/*/metabat/data_bin.4.fa \
    output/data/1/binning/*/metabat/data_bin.5.fa \
    output/data/1/binning/*/metabat/data_bin.6.fa \
    output/data/1/binning/*/metabat/data_bin.7.fa \
    output/data/1/binning/*/metabat/data_bin.8.fa \
    output/data/1/binning/*/metabat/data_bin.9.fa \
    output/data/1/binning/*/metabat/data_bin.10.fa \
    output/data/1/binning/*/metabat/data_notBinned.fa
    ```

!!! Question "Task 6"
    Now inspect the QUAST report with:
    ```BASH
    firefox quast_bins/report.html
    ```
    Can you identify, where large part of the genome fractions are assigned to? What might be the reason? Have a close look on the contig statistics of each genome bin and the unbinned fraction
    ??? Solution
        Large parts of the contigs have been assigned to the unbinned fraction. The reason for that is the small contig size. Metabat2 uses 2500 bp as default cutoff for contigs that can be assigned to genome bins, other contigs are automatically assigned to the unbinned fraction. This cutoff can be lowered to 1500 bp but would not improve the results very much here, since so many contigs are smaller than that.

In summary, three of your bins meet at least the completeness and contamination criteria to be considered as medium quality (Completeness > 50% and Contamination < 10%) 
according to the [**Minimum Information about a Metagenome-Assembled Genome (MIMAG)**](https://www.gensc.org/pages/standards/checklists.html) standard. 
In the next section we will examine their taxonomy.

---

➡️ [**Continue to: Classification**](./classification.md) 
