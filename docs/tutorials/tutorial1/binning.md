We will now run a binning using the mappings results from the previous section

##Metabat

MetaBAT, An Efficient Tool for Accurately Reconstructing Single
Genomes from Complex Microbial Communities.

Grouping large genomic fragments assembled from shotgun metagenomic
sequences to deconvolute complex microbial communities, or metagenome
binning, enables the study of individual organisms and their
interactions. MetaBAT is an automated metagenome binning software
which integrates empirical probabilistic distances of genome abundance
and tetranucleotide frequency. See the `MetaBAT home page
<https://bitbucket.org/berkeleylab/metabat>`_
for more info.
  
We will now run a binnning using metabat, we add these lines below the assembly part of the parameter file (the bwa part has already been added, don't duplicate it):
```BASH
---8<--- "default/tutorials/tutorial1/fullpipeline_binning.yml:39:52"
```
The complete parameter file is appended at the end of this page.

We can run it directly with:
```BASH
---8<--- "scripts/tutorials/tutorial1/test_binning.sh:3:12"
```
  
  
MetaBAT will generate 10 bins from our assembly::
```
  ls -l output/data/1/binning/0.5.0/metabat/
  data_bin.1.fa
  data_bin.10.fa
  data_bin.2.fa
  data_bin.3.fa
  data_bin.4.fa
  data_bin.5.fa
  data_bin.6.fa
  data_bin.7.fa
  data_bin.8.fa
  data_bin.9.fa
  data_bin_contig_mapping.tsv
  data_bins_stats.tsv
  data_contigs_depth.tsv
  data_notBinned.fa
```

You can have a look on some bin statistics:
```
  less output/data/1/binning/0.5.0/metabat/data_bins_stats.tsv
```

=====
The complete parameter file:
```BASH
---8<--- "default/tutorials/tutorial1/fullpipeline_binning.yml"
```




