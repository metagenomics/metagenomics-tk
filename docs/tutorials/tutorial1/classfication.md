Taxonomonic classification tools assign taxonommic labels to reads or assembled contigs of metagenomic datasets.
We will perform taxonomic classification of our genome bins using GTDBtk.

## Genome Taxonomy Database (GTDB)

GTDB-Tk is a software toolkit for assigning objective taxonomic 
classifications to bacterial and archaeal genomes based on the 
[Genome Taxonomy Database GTDB](https://gtdb.ecogenomic.org>). 
It is designed to work with recent advances that allow hundreds or thousands of metagenome-assembled 
genomes (MAGs) to be obtained directly from environmental samples. 
It can also be applied to isolate and single-cell genomes. 

See the [GTDBTk homepage](https://ecogenomics.github.io/GTDBTk/index.html) for more info.


Next, let's assign taxonomic labels to our binning results using GTDB-TK. 

The following lines have to added to our parameter file in order to run the classification:
```BASH
---8<--- "default/tutorials/tutorial1/fullpipeline_classification.yml:60:66"
```
TODO: passt das oben?

!!! Question "Task 1"
    Run the following command to for the classification:
    ```BASH
    ---8<--- "scripts/tutorials/tutorial1/test_classification.sh:3:12"
    ```
TODO: Datenbank Pfad unbedingt in diesem Aufruf anpassen! 


!!! Question "Task 2"
    Why was the classification so fast?
    !!! Solution
        The classification was quite fast due to mash/ANI fast identification (see https://ecogenomics.github.io/GTDBTk/commands/classify_wf.html). Marker gene based classification did not happen, in a real world metagenome sample, the classification would take much longer.
  
!!! Question "Task 3"
    Locate the classification results inside your `output` directory and find out the classfication for each bin.
    ??? Solution
        The results are stored in `TODO: Verzeichnis!!`. You can find the classification results in the file `TODO: file`:
        ```BASH
        cd ~/mgcourse/
        TODO: anpassen
        cut -f?,? ???
        ...
        ```
!!! Question "Task 4"
    Inspect the results in EMGB!

For reference, your complete parameter file should look like this:
??? Parameter-file

    ```BASH
    ---8<--- "default/tutorials/tutorial1/fullpipeline_classification.yml"
    ```      



