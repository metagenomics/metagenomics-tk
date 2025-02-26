To get an impression about the quality of our bins, we compute the completeness and contamination values for our bins. 

## Computing completeness and contamination using CheckM

CheckM provides a set of tools for assessing the quality of genomes recovered from isolates, single cells, or metagenomes. It provides robust estimates of genome completeness and contamination by using collocated sets of genes that are ubiquitous and single-copy within a phylogenetic lineage. Assessment of genome quality can also be examined using plots depicting key genomic characteristics (e.g., GC, coding density) which highlight sequences outside the expected distributions of a typical genome. CheckM also provides tools for identifying genome bins that are likely candidates for merging based on marker set compatibility, similarity in genomic characteristics, and proximity within a reference genome tree.
See the `CheckM home page <https://ecogenomics.github.io/CheckM/>`_ for more info.

We will now CheckM2 using the Metagenomics Toolkit - we add these lines below the binning part of the parameter file:
```BASH
---8<--- "default/tutorials/tutorial1/fullpipeline_bin_quality.yml:53:59"
```
The complete parameter file is appended at the end of this page.

We can run it directly with:
```BASH
---8<--- "scripts/tutorials/tutorial1/test_bin_quality.sh:3:12"
```

Now, the results can be found here::
```
  cd /vol/WGS-data/
  less output/data/1/magAttributes/3.0.0/checkm2/data_checkm2_generated.tsv
```

=====
The complete parameter file:
```BASH
---8<--- "default/tutorials/tutorial1/fullpipeline_bin_quality.yml"
```

