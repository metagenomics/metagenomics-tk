# Meta-Omics-Toolkit

The Meta-Omics-Toolkit allows you to run the full pipeline of assembly, binning and many other downstream analysis
tasks.

All tools follow the same error strategy. The execution of a tool is retried three times. If the run fails the fourth time,
it will be ignored. If the execution is ignored, the toolkit will continue to run all tools that do not depend on the output
of the failed tool run. Exceptions of this handling are specified in the corresponding module section.

## Used databases and tools
For the annotation of metagenome assembled genomes (MAGs), plasmid identification and other MAG the meta-omics-toolkit uses the following databases and tools:
  - GTDB (Genome Taxonomy Database)
  - CheckM
  - VFDB (Virulence Factors Database)
  - KeGG
  - Prokka
  - bacmet20
  - uniref90
  - Comprehensive Antibiotics Resistance Database (CARD)
  - SCAPP for plasmid assembly
  - PLSD - a plasmid database

## Modules
### Quality Control
The quality control module removes adapters, trims and filters short read and long read data.

### Assembly and Binning
The assembly and binning module uses the ultra-fast de novo assembler Megahit(Li et al, 2015) for large and complex metagenomics assembly via succinct de Bruijn graphs and MetaBAT
(Kang et al. 2019) for robust and efficient genome reconstruction from assembled metagenomes.

### Annotation and Mag Attributes
The annotation module is able to predict genes and annotate those based on Prokka and a set of user provided databases.
In addition, the resistance gene identifier is executed by default.

### Plasmids
The plasmid module is able to identify contigs as plasmids and also to assemble plasmids from the samples fastq data.
The module is executed in two parts. In the first part contigs of a metagenome assembler are scanned for plasmids. In the second
part a plasmid assembler is used to assemble circular plasmids out of raw reads. All plasmid detection tools are executed on the
circular assembly result and on the contigs of the metagenome assembler. Just the filtered sequences are used for downstream analysis.

The identification of plasmids is based on the combined result of tools which have a `filter` property assigned. Results of all
tools that have the `filter` property set to true are combined either by a logical `OR` or by a logical `AND`.

Example for the `OR` and `AND` operations: Let's assume that we have three plasmid detection tools (t1, t2, t3) that have four contigs
(c1, c2, c3, c4) as input. Let's further assume that c1 and c2 are detected by all tools as contigs and c3 and c4 are only detected
by t1 and t2. By using an `AND` only c1 and c2 are finally reported by the module as plasmids. By using an `OR` all contigs would be
annotated as plasmids.

Only the detected plasmids will be used for downstream analysis. The read mapper can either be Bowtie or Bwa for Illumina and minimap
for long reads.

## References
  - Li, D., Liu, C-M., Luo, R., Sadakane, K., and Lam, T-W., (2015) MEGAHIT: An ultra-fast single-node solution for large and complex metagenomics assembly via succinct de Bruijn graph. Bioinformatics, doi: [10.1093/bioinformatics/btv033](https://doi.org/10.1093/bioinformatics/btv033 ) [PMID: 25609793].
  - Kang DD, Li F, Kirton E, Thomas A, Egan R, An H, Wang Z. MetaBAT 2: an adaptive binning algorithm for robust and efficient genome reconstruction from metagenome assemblies. PeerJ. 2019 Jul 26;7:e7359. doi: [10.7717/peerj.7359](https://doi.org/10.7717/peerj.7359). PMID: 31388474; PMCID: PMC6662567.
