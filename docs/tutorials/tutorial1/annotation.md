## Annotation

The Metagenomics-Toolkit provides several tools for the annotation of genes, the basis for them is the gene prediction
which is done using prodigal. Prodigal is also part of prokka, which also creates a quick functional annotation of our data. 
Annotations with larger databases are out of the scope of this workshop due to runtime limitations, but we will inspect them
nevertheless within EMGB.

## Prokka

Whole genome annotation is the process of identifying features of interest in a set of genomic DNA sequences, 
and labelling them with useful information. Prokka is a software tool to annotate bacterial, archaeal and viral 
genomes quickly and produce standards-compliant output files. See [Prokka homepage](https://github.com/tseemann/prokka) 
for more info.

The following lines have to added to our parameter file in order to run the annotation with prokka:
```BASH
---8<--- "default/tutorials/tutorial1/fullpipeline_annotation.yml:67:75"
```
TODO: passt das oben?

!!! Question "Task 1"
    Run the following command to for the annotation:
    ```BASH
    ---8<--- "scripts/tutorials/tutorial1/test_annotation.sh:3:12"
    ```
!!! Question Task 2
    Locate the annotation files inside the `output` directory and inspect the `.gff.` file.
    ??? Solution
        The 
        ```BASH
        cd ~/mgcourse
        less ?...gff
        ```
TODO: Pfad anpassen hier oben

!!! Question Task 3
    Next, inspect the annotation results in EMGB.
