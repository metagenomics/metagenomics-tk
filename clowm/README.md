# Metagenomics-Toolkit

[![DOI:10.1101/2024.10.22.619569](http://img.shields.io/badge/DOI-10.1101/2024.10.22.619569-B31B1B.svg)](https://doi.org/10.1101/2024.10.22.619569)

The Metagenomics-Toolkit is a scalable, data agnostic workflow that automates the analysis of short and long metagenomic reads obtained from Illumina or Oxford Nanopore Technology devices, respectively.
The Toolkit offers not only standard features expected in a metagenome workflow, such as quality control, assembly, binning, and annotation, but also distinctive features,
such as plasmid identification based on various tools, the recovery of unassembled microbial community members, and the discovery of microbial interdependencies through a combination of dereplication, co-occurrence, and genome-scale metabolic modeling.
Furthermore, the Metagenomics-Toolkit includes a machine learning-optimized assembly step that tailors the peak RAM value requested by a metagenome assembler to match actual requirements, thereby minimizing the dependency on dedicated high-memory hardware.

**Schema of the complete Metagenomics-Toolkit workflow:**

![per-sample-workflow](https://openstack.cebitec.uni-bielefeld.de:8080/clowm/IBG-5_Grafik-Veröffentlichung-A4_V11.png)

**Below is a list of tools and databases that are enabled for the CloWM service.**
**Only the best practice tools are enabled. More tools and databases will be enabled in the future.**

## Modules and Tools

Currently the following modules and tools are enabled for execution in CloWM:

| Box number in figure | Modules | Tools |
|-----|----|-------|
| 2 | Quality Control | Fastp (Illumina), Porechop (Nanopore), Filtlong (Nanopore), KMC, Nonpareil |
| 3 | Assembly | MEGAHIT (Illumina), Metaflye (Nanopore) |
| 4 | Read Mapping | BWA-MEM2 (Illumina), Minimap (Nanopore) |
| 6 | Binning | MetaBAT2, MAGScoT |
| 9 | Plasmids Assembly /Examination | Platon, ViralVerify, PlasClass, PLSDB, SCAPP |
| 8 | Phylogeny/Taxonomy and Annotation | GTDB-tk, CheckM2, Prokka, RGI, MMseqs2, MMSeqs2 taxonomy |

## Databases

The following databases are used:
  - GTDB (Genome Taxonomy Database)
  - VFDB (Virulence Factors Database)
  - KEGG (Kyoto Encyclopedia of Genes and Genomes)
  - bacmet20 (Antibacterial Biocide- and Metal Resistance Genes Database)
  - uniref90 (UniProt Reference Cluster Database)
  - CARD (Comprehensive Antibiotics Resistance Database)
  - PLSDB (Plasmid database)

## Documentation

Further documentation regarding the CLI tool and module description can be found [here](https://metagenomics.github.io/metagenomics-tk)

## Citations

If you use this pipeline, please cite:

Belmann, P., Osterholz, B., Kleinboelting, N., Puehler, A., Schlueter, A., & Sczyrba, A. (2024). 
Metagenomics-Toolkit: The Flexible and Efficient Cloud-Based Metagenomics Workflow featuring Machine Learning-Enabled Resource Allocation. 
Cold Spring Harbor Laboratory. https://doi.org/10.1101/2024.10.22.619569 
