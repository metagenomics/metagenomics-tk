Annotation is the process of identifying features of interest in a set of genomic DNA sequences and labeling them with information such as their function.
The Metagenomics-Toolkit provides several tools for annotating genes, based on the gene prediction using Prodigal.
Prodigal is also part of Prokka, which provides a fast functional annotation of our data in addition to gene prediction. 

Annotation with larger databases is beyond the scope of this workshop due to runtime limitations, but we will explore it within EMGB.

## Prokka

Prokka is an efficient, user-friendly and open source bioinformatics tool designed for the annotation of bacterial genomes.
It automates the prediction of genes, tRNAs, rRNAs, and other genomic features, utilizing various databases and algorithms to ensure accurate annotations.
Prokka supports standard output formats such as GenBank and GFF, facilitating further analysis with compatible tools. 

See [Prokka homepage](https://github.com/tseemann/prokka) for more information.

The following snippet represents the Toolkit contiguration for the annotion module, just running Prokka:
```YAML linenums="1" title="Annotation Configuration File Snippet 1"
---8<--- "default/tutorials/tutorial1/fullpipeline_annotation.yml:75:78"
```

!!! Question "Task 1"
    Run the following command to for the annotation:
    ```BASH
    ---8<--- "scripts/tutorials/tutorial1/test_annotation.sh:3:13"
    ```

For reference, your complete parameter file should look like this:
??? Parameter-file

    ```BASH
    ---8<--- "default/tutorials/tutorial1/fullpipeline_annotation.yml"
    ```

!!! Question "Task 2"
    Locate the annotation files in the `output` directory and look for `*.txt` files that summarize the number of genes per category.  
    Can you report a bin that contains CRISPR sequences?
    ??? Solution
        ```BASH
        cd ~/mgcourse
        ```
        In the following case a bin with the number 3 contains CRISPR sequences. In your case maybe another bin contains CRISPR sequences.
        ```BASH
        cat output/data/1/annotation/1.0.0/prokka/data_bin.3.txt
        ```
        The above command leads to the following output for the pre computed results: 
        ```BASH
        organism: Genus species strain 
        contigs: 150
        bases: 520975
        CDS: 853
        CRISPR: 3
        tRNA: 12
        ```

If you want to look into the descriptions of the annotated genes then you have to look into the GFF (General Feature Format) files. 

A GFF file is a tab-delimited text format used in bioinformatics to describe genomic features such as genes, exons, introns, and more.
It consists of multiple columns that provide specific details about each feature:

  * Sequence Name: The identifier for the sequence or chromosome.

  * Source: The source or database providing the annotation.

  * Feature Type: The type of feature (e.g., gene, exon, CDS).

  * Start and End Positions: The coordinates where the feature begins and ends.

  * Score: A confidence score for the feature prediction.

  * Strand: Indicates the directionality (forward or reverse) of the feature.

  * Phase: For features like CDS, it indicates reading frame offsets.

  * Attributes: Additional information such as gene IDs or product names.

Optional header lines starting with '##' provide metadata about the file, adding context to the annotations. 

!!! Question "Task 3"
    Find a gff file of a bin that contains tRNAs the same way as you did in the previous taks. Search with `zgrep -P "\ttRNA\t"` for tRNAs and the corresponding annotations file (gff).
    Can you tell the products of two tRNAs. 

    ??? Solution
        In the our precomputed case bin3 contains two rRNAs.
        ```
        zgrep tRNA output/data/1/annotation/1.0.0/prokka/data_bin.3.gff.gz
        ```
        For example, there the tRNA products `tRNA-Thr(cgt)` and `tRNA-His(gtg)`.  
        

### EMGB

Using the Toolkit, you have now predicted and annotated hundreds of genes.
If you were to use all the tools in the annotation module, you could annotate even more genes using databases such as KEGG.
While you can of course use the command line to search for genes,
MAGs, and pathways in this mass of data, you could also explore your dataset with EMGB.

!!! Question "Task 4"
    Lets imagine you search for MAGs that contain genes that encode for the enzyme `5.3.1.9`.
    Hint: The enzyme is involved the glycolysis pathway. 

    Click on the enzyme and a filter will be automatically created. Now if you click on the "MAGs" tab, you should see MAGs that contain that gene.
    How many MAGs contain that gene?
    Which MAG has the highest completeness?

    In the next step, you might also interested in neighboring genes of the MAG with the highest completeness. What is the name of the gene to the left of it?
    
    ??? Solution
        In the precomputed output the results look the following way. In your case it might another bin id might have the highest completeness.
        Five MAGs contain the gene that encodes the enzyme `5.3.1.9`. The MAG with the highest completeness is `Bin 6`. To the left of the gene encoding the enzyme is the gene for `Proline--tRNA ligase`.
