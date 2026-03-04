In this part of the tutorial we evaluate the recovered Bins for quality, classify them taxonomically and annotate their functional content.
First, we run CheckM2, the machine‑learning based version of CheckM that assesses completeness and contamination by leveraging the full repertoire of genomic markers—including multi‑copy genes
and metabolic modules—rather than relying solely on single‑copy core genes. Next, we use GTDB‑Tk to assign each Bin to the Genome Taxonomy Database (GTDB),
a regularly updated, phylogenetic classification of bacteria and archaea that employs ANI and phylogenetic markers to place query genomes within a standardized taxonomy.
We will then annotate the genomes with Prodigal (via Prokka) to predict genes and rapidly assign functions, providing a comprehensive view of each Bin’s biological potential.
Finally, we apply RGI to examine the Bins with regard to antibiotic resistance genes.

Together, these steps deliver a complete assessment of Bin quality, lineage, and gene function within the Metagenomics‑Toolkit workflow.

## Metagenomics-Toolkit Execution

The following lines represent the part of the configuration that tells the Toolkit to run the MagAttributes and Annotation module that include the CheckM2, GTDB-Tk, Prokka and Resistance Gene Identifier (RGI) tool:

```YAML linenums="1" title="MagAttributes Annotation Configuration File Snippet 1"
---8<--- "default/tutorials/tutorial2/fullPipeline_lineage_and_function.yml:69:93"
```
    
For reference, your complete parameter file looks like this:
??? Parameter-file

    ```BASH
    ---8<--- "default/tutorials/tutorial2/fullPipeline_lineage_and_function.yml"
    ```    

If you run the Toolkit again with the `-resume` flag, it will run all additional modules and tools you have specified.

!!! Question "Task 1"

    Change into the directory of the Toolkit session.
    ```BASH
    cd /vol/volume/sessions/metagenomics_metagenomics-tk 
    ```

    Run the MagAttributes annotation module with the Metagenomics-Toolkit using the following command:
    ```BASH
    ---8<--- "scripts/tutorials/tutorial2/test_lineage_and_function.sh:5:14"
    ```
    This command will take between four and six minutes to run.

## Output

### Computing Completeness and Contamination using CheckM2

!!! Question "Task 2"
    Locate the CheckM results inside the `output/SRR492065/` directory and find out the completeness and contamination values for all of our Bins.
    ??? Solution
        Change into the directory of the Toolkit session.
        ```BASH
        cd /vol/volume/sessions/metagenomics_metagenomics-tk 
        ```
        View the Bin ID, Bin completeness and contamination:
        ```BASH
        cut -f2,3,4 output/SRR492065/1/magAttributes/*/checkm2/SRR492065_checkm2_generated.tsv
        ```
        Output:
        ```BASH        
        BIN_ID  COMPLETENESS    CONTAMINATION
        SRR492065_bin.1.fa      54.16   0.31
        SRR492065_bin.2.fa      49.86   0.15
        SRR492065_bin.3.fa      94.75   0.09
        SRR492065_bin.4.fa      85.99   18.95
        SRR492065_bin.5.fa      96.32   0.14
        ``` 
        
According to Bowers et al. [^1], Bins can be grouped into three different quality categories.
Part of this quality assessment involves evaluating Bin completeness and contamination.

* High quality: > 90% complete and < 5% contamination 

* Medium quality: >= 50% complete and < 10% contamination 

* Low quality: < 50% complete and < 10% contamination 

Based on this definition we can conclude that there are 2 high quality, 1 medium quality and 1 low quality Bin.

### Genome Taxonomy Database (GTDB)

Next the Toolkit performed taxonomic classification of our genome Bins using GTDB-Tk.
GTDB-Tk is a software toolkit that enables users to classify their own genomic or metagenomic data against the GTDB taxonomy.

See the [GTDB-Tk homepage](https://ecogenomics.github.io/GTDBTk/index.html) for more information.

The classification itself is quite fast. The fast execution can be explained by the fast identification of a genome, and thus a taxonomic label, in the GTDB via Mash and ANI.
Classification based on marker genes was not necessary. In a real metagenome sample, the classification would usually take longer.

!!! Question "Task 3"
    Check the classification and BIN_ID column of the gtdb output for the SRR492065 sample. What is the classification of genomes that are at least 50% complete and are at most 10% contaminated?
 
    ??? Solution
        Change into the directory of the Toolkit session.

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
        SRR492065_bin.2.fa  d__Bacteria;p__Bacillota;c__Bacilli;o__Lactobacillales;f__Enterococcaceae;g__Enterococcus;s__Enterococcus faecalis
        SRR492065_bin.3.fa  d__Bacteria;p__Actinomycetota;c__Actinomycetes;o__Propionibacteriales;f__Propionibacteriaceae;g__Cutibacterium;s__Cutibacterium avidum
        SRR492065_bin.4.fa  d__Bacteria;p__Bacillota;c__Bacilli;o__Staphylococcales;f__Staphylococcaceae;g__Staphylococcus;s__Staphylococcus aureus
        SRR492065_bin.5.fa  d__Bacteria;p__Bacillota;c__Clostridia;o__Tissierellales;f__Peptoniphilaceae;g__Peptoniphilus_A;s__Peptoniphilus_A lacydonensis
        ```

        ```TSV
        BIN_ID  COMPLETENESS    CONTAMINATION
        SRR492065_bin.1.fa      54.16   0.31
        SRR492065_bin.2.fa      49.86   0.15
        SRR492065_bin.3.fa      94.75   0.09
        SRR492065_bin.4.fa      85.99   18.95
        SRR492065_bin.5.fa      96.32   0.14
        ```

        Based on the output and the comparison of the BIN_ID columns, we can say that the following species could be detected that belong to bins that are at least 50% complete and at most 10% contaminated:
        ```BASH
        d__Bacteria;p__Bacillota;c__Bacilli;o__Lactobacillales;f__Enterococcaceae;g__Enterococcus;s__Enterococcus faecalis
        d__Bacteria;p__Actinomycetota;c__Actinomycetes;o__Propionibacteriales;f__Propionibacteriaceae;g__Cutibacterium;s__Cutibacterium avidum
        d__Bacteria;p__Bacillota;c__Clostridia;o__Tissierellales;f__Peptoniphilaceae;g__Peptoniphilus_A;s__Peptoniphilus_A lacydonensis
        ```

You now have detected two pathogens and since many bacterial pathogens routinely carry antibiotic‑resistance genes we now
will now investigate their genes using Prokka and RGI.

### Annotate Genes with Prokka and RGI

Prokka is an efficient, user-friendly and open source bioinformatics tool designed for the annotation of bacterial genomes.
Prokka supports standard output formats such as GenBank and GFF, facilitating further analysis with compatible tools. 

See [Prokka homepage](https://github.com/tseemann/prokka) for more information.

We used the genes predicted by Prodigal to identify possible antibiotic resistance genes with the [Resistance Gene Identifier (RGI)](https://github.com/arpcard/rgi) tool.
You can inspect the output of RGI by running the following commands: 

!!! Question "Task 4"
    Change into the directory of the Toolkit session.
    ```BASH
    cd /vol/volume/sessions/metagenomics_metagenomics-tk 
    ```

    List all files that were produced for the SRR492065 sample.
    ```BASH
    ls output/SRR492065/1/annotation/2.0.1/rgi/
    ``` 
    
The tool produces among other things `*.csv` files that contain a summary of all hits found according to the CARD database for every bin. 

!!! Question "Task 5"
    Depending on the bin ID of your Enterococcus faecalis, count the number of antibiotic resistance genes that could be detected in that bin of the sample SRR492065.

    ??? Solution
        Change into the directory of the Toolkit session.
        ```BASH
        cd /vol/volume/sessions/metagenomics_metagenomics-tk 
        ```
        In the following case the Enterococcus faecalis genome has the bin number 2.
        ```BASH
        cat output/SRR492065/1/annotation/2.0.1/rgi/SRR492065_bin.2.fa.rgi-1.csv
        ```

        In the output you can count the detected genes (6):
        ```BASH
        gene,SRR492065_bin.2.fa.rgi
        IreK,1
        dfrE,1
        efrA,1
        efrB,1
        lsaA,1
        vanY gene in vanB cluster,1
        ```

---

➡️ [**Continue to:Dereplication**](./dereplication.md)

[^1]: https://www.nature.com/articles/nbt.3893
