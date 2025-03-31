Taxonomic classification tools assign taxonomic labels to reads or assembled contigs of metagenomic datasets.
We will perform taxonomic classification of our genome bins using GTDB-Tk.

## Genome Taxonomy Database (GTDB)

[GTDB](https://gtdb.ecogenomic.org) is an open-access database that provides a comprehensive taxonomy for bacterial and archaeal genomes.
It uses genome-scale phylogenetic analysis to classify microbial genomes into a standardized taxonomic framework.
The database is regularly updated with new genomes, ensuring it reflects the latest advancements in microbial taxonomy.
GTDB-Tk is a software toolkit that enables users to classify their own genomic or metagenomic data against the GTDB taxonomy.
It uses metrics like average nucleotide identity (ANI) and other phylogenetic markers to determine the taxonomic placement of query genomes.
This tool is particularly useful for researchers who want to assign taxonomic labels to their sequences within a consistent and widely-accepted framework,
facilitating better understanding of microbial diversity and evolution.

See the [GTDB-Tk homepage](https://ecogenomics.github.io/GTDBTk/index.html) for more information.

Next, let's assign taxonomic labels to our binning results using GTDB-Tk. 
The following snippet represents the Toolkit configuration for the classification part of the MagAttributes module:
```YAML linenums="1" title="Classification Configuration File Snippet 1"
---8<--- "default/tutorials/tutorial1/fullpipeline_classification.yml:68:74"
```

!!! Question "Task 1"
    Run the following command for the classification:
    ```BASH
    ---8<--- "scripts/tutorials/tutorial1/test_classification.sh:3:13"
    ```

For reference, your complete parameter file should look like this:
??? Parameter-file

    ```BASH
    ---8<--- "default/tutorials/tutorial1/fullpipeline_classification.yml"
    ```

The classification took some time, most of it is due to a mash database, that needs to be built before the actual classification. After that, the classification itself was quite fast. The fast execution can be explained by the fast identification of a genome, and thus a taxonomic label, in the GTDB via Mash and ANI.
Classification based on marker genes was almost not necessary. In a real metagenome sample, the classification would usually take longer.
You can read [here](https://ecogenomics.github.io/GTDBTk/commands/classify_wf.html) about the different steps of the GTDB-Tk classification.

!!! Question "Task 2"
    Inspect the different columns of the GTDB-Tk output file with the following command to see how many genomes were processed by the
    ANI screening and how many were identified by the marker gene approach.
     
    ```
    cd ~/mgcourse/
    ```

    ```
    cat output/data/1/magAttributes/*/gtdb/data_gtdbtk_generated_combined.tsv  | column -s$'\t' -t | less -S
    ```

    ??? Solution
        For example, the column `fastani_ani` reports that 9 out of 10 genomes were processed using ANI. 
        ```TSV
        fastani_ani
        98.52
        97.43
        98.4
        98.7
        99.23
        N/A
        99.51
        98.58
        99.67
        98.82
        ```

Now, of course, you are interested in the classification of your genomes. The results can be affected by missing marker genes or contamination, as indicated on the GTDB-Tk [website](https://ecogenomics.github.io/GTDBTk/commands/classify_wf.html).
You already have an estimate of the completeness and contamination of your genomes from the bin quality sections and can use this information for the next task.
  
!!! Question "Task 3"
    Check the classification and BIN_ID column of the gtdb output. What is the classification of genomes that are at least 50% complete and are at most 10% contaminated?
 
    ??? Solution
        ```BASH
        cd ~/mgcourse/
        ```

        By executing the following two commands you get the classification and the quality values:
        ```BASH
        cut -f 1,5 output/data/1/magAttributes/*/gtdb/data_gtdbtk_generated_combined.tsv | column -s$'\t' -t
        ```

        ```BASH
        cut -f2,3,4 output/data/1/magAttributes/*/checkm2/data_checkm2_generated.tsv | column -s$'\t' -t
        ```

        The output is the following:
        ```TSV
        BIN_ID          classification
        data_bin.1.fa   d__Bacteria;p__Desulfobacterota;c__Desulfovibrionia;o__Desulfovibrionales;f__Desulfovibrionaceae;g__Lawsonia;s__Lawsonia intracellularis
        data_bin.10.fa  d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Porphyromonadaceae;g__Porphyromonas;s__Porphyromonas gingivalis
        data_bin.2.fa   d__Bacteria;p__Bdellovibrionota;c__Bdellovibrionia;o__Bdellovibrionales;f__Bdellovibrionaceae;g__Bdellovibrio;s__Bdellovibrio bacteriovorus
        data_bin.3.fa   d__Bacteria;p__Aquificota;c__Aquificae;o__Aquificales;f__Aquificaceae;g__Aquifex;s__Aquifex aeolicus
        data_bin.4.fa   d__Bacteria;p__Bacillota_A;c__Clostridia;o__Tissierellales;f__Peptoniphilaceae;g__Finegoldia;s__Finegoldia magna_H
        data_bin.5.fa   d__Bacteria;p__Campylobacterota;c__Campylobacteria;o__Campylobacterales;f__Helicobacteraceae;g__Helicobacter;s__
        data_bin.6.fa   d__Bacteria;p__Chlamydiota;c__Chlamydiia;o__Chlamydiales;f__Chlamydiaceae;g__Chlamydophila;s__Chlamydophila pneumoniae
        data_bin.7.fa   d__Bacteria;p__Actinomycetota;c__Actinomycetia;o__Mycobacteriales;f__Mycobacteriaceae;g__Mycobacterium;s__Mycobacterium leprae
        data_bin.8.fa   d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Wigglesworthia;s__Wigglesworthia glossinidia_B
        data_bin.9.fa   d__Bacteria;p__Fusobacteriota;c__Fusobacteriia;o__Fusobacteriales;f__Fusobacteriaceae;g__Fusobacterium;s__Fusobacterium nucleatum
        ```

        ```TSV
        BIN_ID          COMPLETENESS  CONTAMINATION
        data_bin.1.fa   16.88         0.0
        data_bin.10.fa  24.04         0.04
        data_bin.2.fa   12.81         0.12
        data_bin.3.fa   20.57         0.01
        data_bin.4.fa   52.73         0.3
        data_bin.5.fa   11.09         0.0
        data_bin.6.fa   81.91         0.46
        data_bin.7.fa   20.68         0.75
        data_bin.8.fa   85.4          0.28
        data_bin.9.fa   39.66         0.62
        ```

        Based on the output and the comparison of the BIN_ID columns, we can say that the following species could be detected that belong to bins that are at least 50% complete and at most 10% contaminated:
        ```BASH
        d__Bacteria;p__Bacillota_A;c__Clostridia;o__Tissierellales;f__Peptoniphilaceae;g__Finegoldia;s__Finegoldia magna_H
        d__Bacteria;p__Chlamydiota;c__Chlamydiia;o__Chlamydiales;f__Chlamydiaceae;g__Chlamydophila;s__Chlamydophila pneumoniae
        d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Wigglesworthia;s__Wigglesworthia glossinidia_B
        ```

!!! Question "Task 4"
    Inspect the results in EMGB - (1) filter for a specific species in the tree on the left side then inspect the MAG tab. (2) Do the same by filtering for a specific MAG (using the filter symbol near the bin id). 
    What do you observe? 
    
    ??? Solution
        Not all genes are classified the same way as the corresponding MAG. Genes are classified using an LCA with mmseqs taxonomy, MAGs are classified using GTDBtk. If in doubt, the GTDBtk classification is probably more accurate.
    

---

➡️ [**Continue to: Annotation**](./annotation.md)
