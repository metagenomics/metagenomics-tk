## FastQC

FastQC aims to provide a simple way to do some quality control checks on raw sequence data coming from high throughput sequencing pipelines.
It provides a modular set of analyses which you can use to give a quick impression of whether your data has any problems of which you should be aware before doing any further analysis.

The main functions of FastQC are

* Import of data from BAM, SAM or FastQ files (any variant)

* Providing a quick overview to tell you in which areas there may be problems

* Summary graphs and tables to quickly assess your data

* Export of results to an HTML based permanent report

* Offline operation to allow automated generation of reports without running the interactive application

See the FastQC home page for more info.

To run FastQC on our data, simply type:

```BASH
fastqc read1.fq.gz read2.fq.gz      
```
Todo: Link the output of the smaller dataset.

The output of fastqc for a large dataset looks like [this](https://s3.bi.denbi.de/cmg/mgcourses/mg2025/qc/DMC_BGA16_4_N_R1_fastqc.html){:target="_blank"}.

## Metagenomics-Toolkit  

The following command runs the quality control module of the Toolkit on one machine.

```BASH
---8<--- "scripts/tutorials/tutorial1/test_qc.sh:3:12"
```

### Output

#### Fastp 

For quality and adapter trimming the Toolkit uses the software fastp which trims reads from 5’ end to 3’ end using a sliding window. 
If the mean quality of bases inside a window drops below a specific q-score, the remaining of the read will be trimmed.
If a read get too short during this trimming, it will be discared.

The fastp output for your data can be viewed [here](https://s3.bi.denbi.de/cmg/mgcourses/qc/data_report.html){:target="_blank"}.

Todo: hint regarding differences.

[Fastp Report](https://s3.bi.denbi.de/cmg/mgcourses/mg2025/qc/DMC_BGA16_4_N_report.html){:target="_blank"}

#### Nonpareil

Nonpareil helps to get an estimation of what to expect from the dataset. Through a redundancy analysis of the reads, 
you get an estimation of the average genome coverage (Todo: explain what it means) and the sequencing effort necessary to theoretically retrieve all
genomes.

[Nonpareil Curves](https://s3.bi.denbi.de/cmg/mgcourses/mg2025/qc/data_nonpareil_curves.pdf){:target="_blank"}
