Taxonomonic classification tools assign taxonommic labels to reads or
assembled contigs of metagenomic datasets.


Genome Taxonomy Database (GTDB)
===============

GTDB-Tk is a software toolkit for assigning objective taxonomic 
classifications to bacterial and archaeal genomes based on the 
`Genome Database Taxonomy GTDB <https://gtdb.ecogenomic.org>`_. 
It is designed to work with recent 
advances that allow hundreds or thousands of metagenome-assembled 
genomes (MAGs) to be obtained directly from environmental samples. 
It can also be applied to isolate and single-cell genomes. 

See the `GTDBTk homepage <https://ecogenomics.github.io/GTDBTk/index.html>`_ 
for more info.

Next, let's assign taxonomic labels to our binning results using GTDB-TK. 
We will now run a binnning using metabat, we add these lines below the assembly part of the parameter file (the bwa part has already been added, don't duplicate it):
```BASH
---8<--- "default/tutorials/tutorial1/fullpipeline_classification.yml:60:66"
```

The complete parameter file is appended at the end of this page.

We can run it directly with:
```BASH
---8<--- "scripts/tutorials/tutorial1/test_classification.sh:3:12"
```

TODO: do something with the results?

TODO: check if the stuff below is still needed --> remove

However, the classification was quite fast due to mash/ANI fast identification (see https://ecogenomics.github.io/GTDBTk/commands/classify_wf.html). Marker gene based classification did not happen, we will repeat the classification with another - non artifical - sample. 
Create a new folder and download the data::

  mkdir /mnt/mgexample/
  cd /mnt/mgexample/
  wget https://openstack.cebitec.uni-bielefeld.de:8080/swift/v1/mgexample/examplesample.tar.gz
  tar -xzvf examplesample.tar.gz

We will do the process for the maxbin results again, step by step::

  gtdbtk classify_wf --extension fasta --cpus 28 --genome_dir examplesample --out_dir gtdbtk_out --mash_db /mnt/release207_v2/mash.msh


Load gtdbtk.backbone.bac120.classify.tree, gtdbtk.bac120.classify.tree.3.tree into the ncbi taxonomy viewer:

https://www.ncbi.nlm.nih.gov/tools/treeviewer/

Convert to itol format::

  gtdbtk convert_to_itol --input_tree gtdbtk_out/classify/gtdbtk.backbone.bac120.classify.tree --output_tree test.itol

Load into:

https://itol.embl.de/upload.cgi


## The complete parameter file
```BASH
---8<--- "default/tutorials/tutorial1/fullpipeline_classification.yml:60:66"
```




