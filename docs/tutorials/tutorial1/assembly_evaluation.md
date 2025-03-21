We will then compare the results of our assembly with those of other assemblers and inspect mappings of reads to the assembly.

##MetaQUAST
QUAST stands for QUality ASsessment Tool. The tool evaluates genome
assemblies by computing various metrics.  You can find all project
news and the latest version of the tool at [sourceforge](http://sourceforge.net/projects/quast).  QUAST utilizes MUMmer,
GeneMarkS, GeneMark-ES, GlimmerHMM, and GAGE. In addition, MetaQUAST
uses MetaGeneMark, Krona tools, BLAST, and SILVA 16S rRNA
database. See the [metaQuast home page](http://quast.sourceforge.net/metaquast/)
for more info.

!!! question "Task 1"
    Copy the pre-computed assembly results to your local directory.
    ```BASH
    cd ~/mgcourse
    ```
    ```BASH
    wget https://s3.bi.denbi.de/cmg/mgcourses/mg2025/assembly_results.tar.gz
    tar -xzvf assembly_results.tar.gz
    ```

!!! question "Task 2"
    Then copy your assembly result from the toolkit run to the assembly_results directory:
    ```BASH
    cd ~/mgcourse/
    ```
    ```BASH
    cp output/data/1/assembly/1.2.1/megahit/data_contigs.fa.gz assembly_results/megahit_out/final.contigs.fa.gz
    gunzip assembly_results/megahit_out/final.contigs.fa.gz
    ```
    
!!! question "Task 3"
    In addition, we need to download some references in order to compare them to our assemblies:
    ```BASH
    cd ~/mgcourse/
    ```
    ```BASH
    wget https://s3.bi.denbi.de/cmg/mgcourses/mg2025/genomes.tar.gz
    tar -xzvf genomes.tar.gz
    ```

To call `metaquast.py` we have to provide reference genomes which
are used to calculate a number of different metrics for evaluation of
the assembly. In real-world metagenomics, these references are usually
not available, of course.

!!! question "Task 4"
    Run metaquast:
    ```BASH
    cd ~/mgcourse/
    ```
    ```BASH
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
    -o quast \
    -l MegaHit,metaSPAdes,Ray_31,Ray_51,velvet_31,velvet_51,idba_ud \
    assembly_results/megahit_out/final.contigs.fa \
    assembly_results/metaspades_out/contigs.fasta \
    assembly_results/ray_31/Contigs.fasta \
    assembly_results/ray_51/Contigs.fasta \
    assembly_results/velvet_31/contigs.fa \
    assembly_results/velvet_51/contigs.fa \
    assembly_results/idba_ud_out/contig.fa
    ```

We will now inspect the QUAST reports.

!!! question "Task 5"
    QUAST generates HTML reports including a number of interactive graphics. To access these reports, locate the html report in the `quast` directory load the reports in your web browser.
    ??? Solution
        ```
        firefox quast/report.html
        ```
!!! question "Task 6"
    Which of the assemblers performs best in terms of coverage of genome fraction of the reference genomes? Which assembly result would you prefer from all of them?
    ??? Solution
        MegaHit, metaSPAdes and idba_ud perform best in terms of covered genome fraction. MegaHit tends to generate larger contigs along with misassemblies, metaSPAdes creates shorter
        contigs with less misassemblies, idba_is quite similar to MegaHit.




