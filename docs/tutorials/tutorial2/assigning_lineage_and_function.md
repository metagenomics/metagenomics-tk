In this part of the tutorial we evaluate the recovered Bins for quality, classify them taxonomically and annotate their functional content.
First, we run CheckM 2, the machine‑learning based version of CheckM that assesses completeness and contamination by leveraging the full repertoire of genomic markers—including multi‑copy genes
and metabolic modules—rather than relying solely on single‑copy core genes. Next, we use GTDB‑Tk to assign each Bin to the Genome Taxonomy Database (GTDB),
a regularly updated, genome‑scale phylogenetic classification of bacteria and archaea that employs ANI and other phylogenetic markers to place query genomes within a standardized taxonomy.
Finally, we annotate the genomes with Prodigal (via Prokka) to predict genes and rapidly assign functions, providing a comprehensive view of each bin’s biological potential.

Together, these steps deliver a complete assessment of Bin quality, lineage, and gene function within the Metagenomics‑Toolkit workflow.

## Metagenomics-Toolkit Execution

The following lines represent the part of the configuration that tells the Toolkit to run the MagAttributes module that includes
the CheckM2, GTDB-Tk, prokka and Resistance Gene Identifier (RGI) tool:

```YAML linenums="1" title="MagAttributes Annotation Configuration File Snippet 1"
---8<--- "default/tutorials/tutorial2/fullPipeline_lineage_and_function.yml:75:99"
```
    
For reference, your complete parameter file looks like this:
??? Parameter-file

    ```BASH
    ---8<--- "default/tutorials/tutorial1/fullPipeline_lineage_and_function.yml"
    ```    

If you run the Toolkit again with the `-resume` flag, it will run all additional modules and tools you have specified.

!!! Question "Task 1"
    Run CheckM2 with the Metagenomics-Toolkit using the following command:
    ```BASH
    ---8<--- "scripts/tutorials/tutorial2/test_lineage_and_function.sh:5:14"
    ```

## Output

### Computing Completeness and Contamination using CheckM2


!!! Question "Task 2"
    Locate the CheckM results inside the `output/SRR492065/` directory and find out the completeness and contamination values for all of our bins.
    ??? Solution
        Change to the directory of the Toolkit session.
        ```BASH
        cd /vol/volume/sessions/metagenomics_metagenomics-tk 
        ```
        View the bin ID, bin completeness and contamination:
        ```BASH
        cut -f2,3,4 output/data/1/magAttributes/*/checkm2/data_checkm2_generated.tsv
        ```
        Output:
        ```BASH        
        BIN_ID  COMPLETENESS    CONTAMINATION
        SRR492065_bin.1.fa      49.86   0.15
        SRR492065_bin.2.fa      94.64   0.11
        SRR492065_bin.3.fa      96.32   0.14
        SRR492065_bin.4.fa      54.16   0.31
        SRR492065_bin.5.fa      85.99   18.95
        ``` 


Explain MIMAG high, medium ...
Taxonomic classification tools assign taxonomic labels to reads or assembled contigs of metagenomic datasets.

### Genome Taxonomy Database (GTDB)


Next the Toolkit performed taxonomic classification of our genome bins using GTDB-Tk.
GTDB-Tk is a software toolkit that enables users to classify their own genomic or metagenomic data against the GTDB taxonomy.

See the [GTDB-Tk homepage](https://ecogenomics.github.io/GTDBTk/index.html) for more information.

The classification itself is quite fast. The fast execution can be explained by the fast identification of a genome, and thus a taxonomic label, in the GTDB via Mash and ANI.
Classification based on marker genes was not necessary. In a real metagenome sample, the classification would usually take longer.

!!! Question "Task 3"
    Check the classification and BIN_ID column of the gtdb output for the SRR492065 sample. What is the classification of genomes that are at least 50% complete and are at most 10% contaminated?
 
    ??? Solution
        Change to the directory of the Toolkit session.

        ```BASH
        cd /vol/volume/sessions/metagenomics_metagenomics-tk 
        ```

        By executing the following two commands you get the classification and the quality values:
        ```BASH
        cut -f 1,5 output/SRR492065/1/magAttributes/*/gtdb/SRR492065_gtdbtk_generated_combined.tsv | column -s$'\t' -t
        ```

        To make a comparison with the CheckM values, list all the completeness and contamination values again.
        ```BASH
        cut -f2,3,4 output/SRR492065/1/magAttributes/*/checkm2/SRR492065_checkm2_generated.tsv | column -s$'\t' -t
        ```

        The output for CheckM and GTDB-Tk is the following:
        ```TSV
        BIN_ID              classification
        SRR492065_bin.1.fa  d__Bacteria;p__Bacillota;c__Bacilli;o__Lactobacillales;f__Enterococcaceae;g__Enterococcus;s__Enterococcus faecalis
        SRR492065_bin.2.fa  d__Bacteria;p__Actinomycetota;c__Actinomycetes;o__Propionibacteriales;f__Propionibacteriaceae;g__Cutibacterium;s__Cutibacterium avidum
        SRR492065_bin.3.fa  d__Bacteria;p__Bacillota;c__Clostridia;o__Tissierellales;f__Peptoniphilaceae;g__Peptoniphilus_A;s__Peptoniphilus_A lacydonensis
        SRR492065_bin.4.fa  d__Bacteria;p__Bacillota;c__Bacilli;o__Lactobacillales;f__Enterococcaceae;g__Enterococcus;s__Enterococcus faecalis
        SRR492065_bin.5.fa  d__Bacteria;p__Bacillota;c__Bacilli;o__Staphylococcales;f__Staphylococcaceae;g__Staphylococcus;s__Staphylococcus aureus
        ```

        ```TSV
        BIN_ID  COMPLETENESS    CONTAMINATION
        SRR492065_bin.1.fa      49.86   0.15
        SRR492065_bin.2.fa      94.64   0.11
        SRR492065_bin.3.fa      96.32   0.14
        SRR492065_bin.4.fa      54.16   0.31
        SRR492065_bin.5.fa      85.99   18.95
        ```

        Based on the output and the comparison of the BIN_ID columns, we can say that the following species could be detected that belong to bins that are at least 50% complete and at most 10% contaminated:
        ```BASH
        d__Bacteria;p__Actinomycetota;c__Actinomycetes;o__Propionibacteriales;f__Propionibacteriaceae;g__Cutibacterium;s__Cutibacterium avidum
        d__Bacteria;p__Bacillota;c__Clostridia;o__Tissierellales;f__Peptoniphilaceae;g__Peptoniphilus_A;s__Peptoniphilus_A lacydonensis
        d__Bacteria;p__Bacillota;c__Bacilli;o__Lactobacillales;f__Enterococcaceae;g__Enterococcus;s__Enterococcus faecalis
        ```


### Annotate Genes with Prokka and RGI

Prokka is an efficient, user-friendly and open source bioinformatics tool designed for the annotation of bacterial genomes.
Prokka supports standard output formats such as GenBank and GFF, facilitating further analysis with compatible tools. 

See [Prokka homepage](https://github.com/tseemann/prokka) for more information.

We used the genes predicted by Prodigal to identify possible antibiotic resistance genes with the [Resistance Gene Identifier (RGI)](https://github.com/arpcard/rgi) tool.
You can inspect the output of RGI by running the following command: 

!!! Question "Task 4"
    Change to the directory of the Toolkit session.
    ```BASH
    cd /vol/volume/sessions/metagenomics_metagenomics-tk 
    ```

    List all files that were produced for the SRR492065 sample.
    ```BASH
    ls output/SRR492065/1/annotation/2.0.1/rgi/
    ``` 
    
The tool produces among other things `*.csv` files that contain a summary of all hits found according to the CARD database for every bin. 

!!! Question "Task 5"
    Depending on the bin ID of your Staphylococcus aureus, count the number of antibiotic resistance genes that could be detected in that bin of the sample SRR492065.

    ??? Solution
        Change to the directory of the Toolkit session.
        ```BASH
        cd /vol/volume/sessions/metagenomics_metagenomics-tk 
        ```
        In the following case the Staphylococcus Aureus genome has the bin number 5.
        ```BASH
        cat output/SRR492065/1/annotation/2.0.1/rgi/SRR492065_bin.5.fa.rgi-1.csv
        ```

        In the output you can count the detected genes (18):
        ```BASH
        gene,SRR492065_bin.5.fa.rgi
        Staphylococcus aureus FosB,1
        Staphylococcus aureus LmrS,1
        Staphylococcus aureus norA,2
        arlR,2
        arlS,2
        kdpD,1
        mepA,1
        mepR,1
        mgrA,2
        norC,1
        sdrM,1
        sepA,1
        tet(38),1
        vanT gene in vanG cluster,1
        ```

---
