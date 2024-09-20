# Metagenomics-Toolkit

The Metagenomics-Toolkit is a scalable, data agnostic workflow that automates the analysis of short and long metagenomic reads obtained from Illumina or Oxford Nanopore Technology devices, respectively.
The Toolkit offers not only standard features expected in a metagenome workflow, such as quality control, assembly, binning, and annotation, but also distinctive features,
such as plasmid identification based on various tools, the recovery of unassembled microbial community members, and the discovery of microbial interdependencies through a combination of dereplication, co-occurrence, and genome-scale metabolic modeling.
Furthermore, the Metagenomics-Toolkit includes a machine learning-optimized assembly step that tailors the peak RAM value requested by a metagenome assembler to match actual requirements, thereby minimizing the dependency on dedicated high-memory hardware.

**Schema of the complete Metagenomics-Toolkit workflow:**

![per-sample-workflow](https://openstack.cebitec.uni-bielefeld.de:8080/swift/v1/clowm/IBG-5_Grafik-Veröffentlichung-A4_V09.jpg)


> [!IMPORTANT]
> Below is a list of tools and databases that are enabled for the CloWM service. 
> Currently only the processing of short reads and the corresponding best practice tools are enabled. More tools and databases will be enabled in the future.

## Modules and Tools

Currently the following modules and tools are enabled for execution in CloWM:

| Box number in figure | Modules | Tools |
|-----|----|-------|
| 2 | Quality Control | Fastp, KMC, Nonpareil |
| 3 | Assembly | MEGAHIT |
| 4 | Read Mapping | BWA-MEM2 |
| 6 | Binning | MetaBAT2, MAGScoT |
| 9 | Plasmids Assembly /Examination |  Platon, ViralVerify, PlasClass, PLSDB, SCAPP |
|  8 | Phylogeny/Taxonomy and Annotation | GTDB-tk, CheckM, Prokka, RGI, MMseqs2, MMSeqs2 taxonomy |

### Plasmids 

The plasmid module is able to identify contigs as plasmids and also to assemble plasmids from the sample's FASTQ data. The module is executed in two parts. In the first part contigs of a metagenome assembler are scanned for plasmids. In the second part a plasmid assembler is used to assemble circular plasmids out of raw reads. All plasmid detection tools are executed on the circular assembly result and on the contigs of the metagenome assembler. Only the filtered sequences are used for downstream analysis.

The identification of plasmids is based on the combined result of tools which have a filter property assigned. The results of all tools with the filter property set to true are combined using either a logical OR or logical AND.

Example of the OR and AND operations: Let's assume that we have three plasmid detection tools (t1, t2, t3) that have four contigs (c1, c2, c3, c4) as input. Let's further assume that c1 and c2 are detected by all tools as contigs and c3 and c4 are only detected by t1 and t2. By using an AND only c1 and c2 are finally reported by the module as plasmids. By using an OR all contigs would be annotated as plasmids.

Only the detected plasmids will be used for downstream analysis.

## Databases

The following databases are used:

For the annotation of metagenome assembled genomes (MAGs), plasmid identification and other MAGs the Metagenomics-Toolkit uses the following databases and tools:
  - GTDB (Genome Taxonomy Database)
  - VFDB (Virulence Factors Database)
  - KEGG (Kyoto Encyclopedia of Genes and Genomes)
  - bacmet20 (Antibacterial Biocide- and Metal Resistance Genes Database)
  - uniref90 (UniProt Reference Cluster Database)
  - CARD (Comprehensive Antibiotics Resistance Database)
  - PLSD (A plasmid database)
