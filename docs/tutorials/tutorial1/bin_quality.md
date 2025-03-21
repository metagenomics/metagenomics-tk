To get an impression about the quality of our bins, we compute the completeness and contamination values for our bins. 

## Computing completeness and contamination using CheckM

CheckM provides a set of tools for assessing the quality of genomes recovered from isolates, single cells, or metagenomes. It provides robust estimates of genome completeness and contamination by using collocated sets of genes that are ubiquitous and single-copy within a phylogenetic lineage. Assessment of genome quality can also be examined using plots depicting key genomic characteristics (e.g., GC, coding density) which highlight sequences outside the expected distributions of a typical genome. CheckM also provides tools for identifying genome bins that are likely candidates for merging based on marker set compatibility, similarity in genomic characteristics, and proximity within a reference genome tree.
See the `CheckM home page <https://ecogenomics.github.io/CheckM/>`_ for more info.

!!! Question "Task 1" 
    We will now run CheckM2 using the Metagenomics Toolkit - we need to add these lines below the binning part of the parameter file:
    ```BASH
    ---8<--- "default/tutorials/tutorial1/fullpipeline_bin_quality.yml:53:59"
    ```
TODO: passt das oben?

The complete parameter file is appended at the end of this page.

!!! Question "Task 2"
    Run CheckM2 directly with:
    ```BASH
    ---8<--- "scripts/tutorials/tutorial1/test_bin_quality.sh:3:12"
    ```
TODO: muss die DB vorher runtergeladen werden, oder ist die klein genug?

!!! Question "Task 3"
    Locate the checkm results inside the `output` directory and find out the completeness and contamination values for all of our bins.
    ??? Solution
        ```
        cd ~/mgcourse
        cut -f ?,? output/data/1/magAttributes/3.0.0/checkm2/data_checkm2_generated.tsv
        TODO: put the results from the table here!!
        ```        

!!! Question "Task 4"
    Compare the CheckM2 results with those from QUAST - do they match your expectations?
    ??? Solution
        No, completeness and contamination are much lower than expected from the covered genome fraction we observed in the QUAST results (see Assembly evaluation part).

We will now run QUAST again on one of the BINs and on the unbinned contigs.

!!! Question "Task 5"
    TODO: Check and test this command!
    Run metaquast with individual BINS as input:
    ```BASH
      cd ~/mgcourse/
   
      metaquast.py --threads 28 --gene-finding \
      -R genomes/Aquifex_aeolicus_VF5.fna,\
      sgenomes/Bdellovibrio_bacteriovorus_HD100.fna,\
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
      output/data/1/binning/0.5.0/metabat/data_bin.1.fa \
      output/data/1/binning/0.5.0/metabat/data_bin.2.fa \
      output/data/1/binning/0.5.0/metabat/data_bin.3.fa \
      output/data/1/binning/0.5.0/metabat/data_bin.4.fa \ 
      output/data/1/binning/0.5.0/metabat/data_bin.5.fa \
      output/data/1/binning/0.5.0/metabat/data_bin.6.fa \
      output/data/1/binning/0.5.0/metabat/data_bin.7.fa \
      output/data/1/binning/0.5.0/metabat/data_bin.8.fa \
      output/data/1/binning/0.5.0/metabat/data_bin.9.fa \
      output/data/1/binning/0.5.0/metabat/data_bin.10.fa \
      output/data/1/binning/0.5.0/metabat/data_notBinned.fa \
    ```

!!! Question "Task 6"
    Now inspect the QUAST report with:
    ```BASH
    firefox quast_bins/report.html
    ```
    Can you identify, where large part of the genome fractions are assigned to? What might be the reason? Have a close look on the contig statistics of each genome bin and the unbinned fraction
    ??? Solution
        Large parts of the the contigs have been assigned to the unbinned fraction. The reason for that is the small contig size, metabat2 uses 2500 bp as default cutoff for contigs that can be assigned to genome bins, other contigs are automatically assigned to the unbinned fraction. This cutoff can be lowered to 1500 bp but would not improve the results very much here, since so many contigs are smaller than that.
TODO: checken, ob man das wirklich in den results sieht!


We will now inspect the QUAST reports.

!!! question "Task 5"
    QUAST generates HTML reports including a number of interactive graphics. To access these reports, locate the html report in the `quast` directory load the reports in your web browser.
    !!! Solution
        ```
        firefox quast/report.html
        ```

For reference, your complete parameter file should look like this:
??? Parameter-file
    ```BASH
    ---8<--- "default/tutorials/tutorial1/fullpipeline_bin_quality.yml"
    ```
