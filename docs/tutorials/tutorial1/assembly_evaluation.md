We will now compare the results of our assembly with those of other assemblers and inspect mappings of reads to the corresponding assemblies.

## MetaQUAST

QUAST stands for QUality ASsessment Tool. The tool evaluates genome
assemblies by computing various metrics.  You can find all project
news and the latest version of the tool at [SourceForge](http://sourceforge.net/projects/quast). 
QUAST is a comprehensive quality assessment tool for genome assemblies that employs several bioinformatics utilities to evaluate assembly accuracy and quality.
It utilizes MUMmer (for genome alignment), GeneMarkS (self-training gene prediction in prokaryotes), GeneMark-ES (gene prediction in eukaryotes),
GlimmerHMM (gene finding using Hidden Markov Models), and GAGE (a pipeline for comparing gene predictions across multiple annotations).

Additionally, MetaQUAST extends QUAST's capabilities to metagenomic assemblies by incorporating tools such as MetaGeneMark (specialized for gene prediction in metagenomes),
Krona Tools (for taxonomic classification visualization), BLAST (for sequence similarity searches against reference databases),
and the SILVA 16S rRNA database (for identifying microbial communities).

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
    Then copy your assembly results from the toolkit run to the assembly_results directory:
    ```BASH
    cd ~/mgcourse/
    ```
    ```BASH
    cp output/data/1/assembly/1.2.1/megahit/data_contigs.fa.gz assembly_results/megahit_out/final.contigs.fa.gz
    gunzip -fd assembly_results/megahit_out/final.contigs.fa.gz
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

To call the `metaquast.py` script, we have to provide reference genomes which are used to calculate a variety of different metrics for evaluating the assembly.
In real-world metagenomics, these references are not available, of course.

!!! question "Task 4"
    Run metaquast:
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
    QUAST generates HTML reports including a number of interactive graphics. To access these reports, locate the html report in the `quast` directory and load the reports in your web browser
    via the following command:
    ```
    firefox quast/report.html
    ```

!!! question "Task 6"
    Which of the assemblers performs best in terms of coverage of genome fraction of the reference genomes? Which assembly result would you prefer among them?
    ??? Solution
        MEGAHIT, metaSPAdes and idba_ud perform best in terms of covered genome fraction. MEGAHIT tends to generate larger contigs along with misassemblies, while metaSPAdes creates shorter
        contigs with fewer misassemblies, idba_ud is quite similar to metaSPades.

---

➡️ [**Continue to: Binning**](./binning.md)
