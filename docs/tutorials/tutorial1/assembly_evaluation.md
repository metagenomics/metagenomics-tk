We will then compare the results of our assembly with those of other assemblers and inspect mappings of reads to the assembly.

##MetaQUAST
QUAST stands for QUality ASsessment Tool. The tool evaluates genome
assemblies by computing various metrics.  You can find all project
news and the latest version of the tool at `sourceforge
<http://sourceforge.net/projects/quast>`_.  QUAST utilizes MUMmer,
GeneMarkS, GeneMark-ES, GlimmerHMM, and GAGE. In addition, MetaQUAST
uses MetaGeneMark, Krona tools, BLAST, and SILVA 16S rRNA
database. See the `metaQuast home page <http://quast.sourceforge.net/metaquast//>`_
for more info.

Copy the pre-computed assembly results to your local directory.
```BASH
  cd /vol/WGS-data
  wget https://s3.bi.denbi.de/cmg/mgcourses/mg2025/assembly_results.tar.gz
  tar -xzvf assembly_results.tar.gz
```
Then copy your assembly result from the toolkit run to the assembly_results directory:
```BASH
  cd /vol/WGS-data/
  cp output/data/1/assembly/1.2.1/megahit/data_contigs.fa.gz assembly_results/megahit_out/final.contigs.fa.gz
  gunzip assembly_results/megahit_out/final.contigs.fa.gz
```

In addition, we need to download some references in order to compare them to our assemblies:
```BASH
  cd /vol/WGS-data/
  wget https://s3.bi.denbi.de/cmg/mgcourses/mg2025/genomes.tar.gz
  tar -xzvf genomes.tar.gz
```

To call ``metaquast.py`` we have to provide reference genomes which
are used to calculate a number of different metrics for evaluation of
the assembly. In real-world metagenomics, these references are usually
not available, of course::
```BASH
  cd /vol/WGS-data
  
  metaquast.py --threads 28 --gene-finding \
  -R /vol/WGS-data/genomes/Aquifex_aeolicus_VF5.fna,\
  /vol/WGS-data/genomes/Bdellovibrio_bacteriovorus_HD100.fna,\
  /vol/WGS-data/genomes/Chlamydia_psittaci_MN.fna,\
  /vol/WGS-data/genomes/Chlamydophila_pneumoniae_CWL029.fna,\
  /vol/WGS-data/genomes/Chlamydophila_pneumoniae_J138.fna,\
  /vol/WGS-data/genomes/Chlamydophila_pneumoniae_LPCoLN.fna,\
  /vol/WGS-data/genomes/Chlamydophila_pneumoniae_TW_183.fna,\
  /vol/WGS-data/genomes/Chlamydophila_psittaci_C19_98.fna,\
  /vol/WGS-data/genomes/Finegoldia_magna_ATCC_29328.fna,\
  /vol/WGS-data/genomes/Fusobacterium_nucleatum_ATCC_25586.fna,\
  /vol/WGS-data/genomes/Helicobacter_pylori_26695.fna,\
  /vol/WGS-data/genomes/Lawsonia_intracellularis_PHE_MN1_00.fna,\
  /vol/WGS-data/genomes/Mycobacterium_leprae_TN.fna,\
  /vol/WGS-data/genomes/Porphyromonas_gingivalis_W83.fna,\
  /vol/WGS-data/genomes/Wigglesworthia_glossinidia.fna \
  -o quast \
  -l MegaHit,metaSPAdes,Ray_31,Ray_51,velvet_31,velvet_51,idba_ud \
  assembly_results/megahit_out/final.contigs.fa \
  assembly_results/metaspades_out/contigs.fasta \
  assembly_results/ray_31/Contigs.fasta \
  assembly_results/ray_51/Contigs.fasta \
  assembly_results/velvet_31/contigs.fa \
  assembly_results/velvet_51/contigs.fa \
  assembly_results/idba_ud_out/contig.fa
```

QUAST generates HTML reports including a number of interactive graphics. To access these reports,
load the reports in your web browser::
```
  firefox quast/report.html
```

## Read Mapping

In this part of the tutorial we will look at the assemblies by mapping
the reads to the assembled contigs.  Different tools exists for
mapping reads to genomic sequences such as `bowtie
<http://bowtie-bio.sourceforge.net/bowtie2/index.shtml>`_ or `bwa
<http://bio-bwa.sourceforge.net/>`_. 

We will now run a mapping using parts of the binning step from the Metagenomics Toolkit, we add these lines below the assembly part of the parameter file:
```BASH
---8<--- "default/tutorials/tutorial1/fullpipeline_map_to_assembly.yml:39:44"
```
The complete parameter file is appended at the end of this page.

We can run it directly with:
```BASH
---8<--- "scripts/tutorials/tutorial1/test_map_to_assembly.sh:3:12"
```

We have to sort the BAM file by starting position of the alignments. This can be done using samtools again::
```BASH
  samtools sort -o output/data/1/binning/0.5.0/contigMapping/data_sorted.bam -@ 28 output/data/1/binning/0.5.0/contigMapping/data.bam 
```
  
Now we have to index the sorted BAM file::
```BASH
  samtools index output/data/1/binning/0.5.0/contigMapping/data_sorted.bam
```
To look at the BAM file use::
```BASH
  samtools view output/data/1/binning/0.5.0/contigMapping/data_sorted.bam | less
```
Now copy/link everything you need for igv in a seperate folder"
```
  cd /vol/WGS-data/
  mkdir igv_data
  cd igv_data
  cp ../output/data/1/assembly/1.2.1/megahit/data_contigs.fa.gz .
  gunzip data_contigs.fa.gz
  ln -s ../output/data/1/binning/0.5.0/contigMapping/data_sorted.bam
  ln -s ../output/data/1/binning/0.5.0/contigMapping/data_sorted.bam.bai
```

We will use the IGV genome browser to look at the mappings::
```BASH
  igv
```


Now let's look at the mapped reads:

1. Load the contig sequences into IGV. Use the menu ``Genomes->Load Genome from File...`` 
2. Load the BAM file into IGV. Use menu ``File->Load from File...`` 

=====
The complete parameter file:
```BASH
---8<--- "default/tutorials/tutorial1/fullpipeline_map_to_assembly.yml"
```




