The tutorial is a concise walkthrough of the analysis workflow that the Metagenomics‑Toolkit provides for metagenomic short‑read projects.
The first step is to run a thorough quality‑control pipeline that trims adapters and low‑quality ends, discards reads that become too short, and removes any host DNA contamination. 
A Nonpareil analysis is then executed to estimate the sequencing depth needed to capture the full microbial diversity.

With the cleaned reads in hand, the Toolkit proceeds to assembly. The MEGAHIT assembler, integrated into the Toolkit, is invoked to build contigs that reflect the underlying genomes present in the sample.
The final part of the analysis is binning. By examining each contig’s k‑mer composition and its coverage, the binning tools cluster contigs that share similar signatures.
These clusters, or “bins,” correspond to draft genomes reconstructed from the community, enabling researchers to recover and study individual organisms without prior knowledge of the sample’s composition.

In short, the guide describes how the Toolkit orchestrates QC, assembly and binning to turn raw metagenomic reads into biologically meaningful genome bins

We will inspect two SRA samples ([SRR492065](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR492065&display=metadata){:target="_blank"} and [SRR492183](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR492183&display=metadata){:target="_blank"}) that belong to the same study "Preborn infant gut metagenome".

## Metagenomics-Toolkit Execution 

The Metagenomics-Toolkit offers multiple tools for each of the aforementioned processing steps.
You will now execute the following Toolkit configuration:

```YAML linenums="1" title="QC,Assembly,Binning Configuration File Snippet 1"
---8<--- "default/tutorials/tutorial2/fullPipeline_reads_to_genomes.yml"
```

!!! question "Task 1"

    You saw the introduction part of the tutorial regarding the structure of the config file. 
    Can you tell which tools or methods are executed in the QC part?

    ??? Solution 
        Fastp, nonpareil, filterHuman, kmc
    
The following configuration runs the tools 

!!! question "Task 2"

    Copy the following command to execute the Toolkit. The Toolkit will need about 3 to 5 (Todo: update) minutes to complete.

    ```BASH
    ---8<--- "scripts/tutorials/tutorial2/test_reads_to_genomes.sh:3:13"
    ```


## Output

In the following we will inspect the outputs of each analysis step.

### Quality Control 

#### Fastp 

Fastp is an efficient, versatile, open-source tool for preprocessing FASTQ files and quality control in next-generation sequencing (NGS) workflows.
Written in C++, fastp provides high-speed performance without compromising accuracy, making it suitable for handling large sequencing data sets.
For quality and adapter trimming fastp trims reads from 5’ end to 3’ end using a sliding window. 

If the mean quality of bases inside a window drops below a specific q-score, the remainder of the read will be trimmed.
If a read gets too short during this trimming, it will be discarded.

!!! info "Quality Report"
    The Metagenomics-Toolkit allows you to run fastp by only reporting the quality of the data. This way you can easily inspect the data
    before actually further processing it. You can do this by setting the `reportOnly` parameter to `true` in the config below. 


!!! question "Task 3"

    You can view the fastp output in the qc output directory:
    ```BASH
    ls -1 output/*/1/qc/*/fastp/
    ```

    which results in the following output:

    ```
    output/SRR492065/1/qc/0.4.0/fastp:
    SRR492065_fastp.json
    SRR492065_fastp_summary_after.tsv
    SRR492065_fastp_summary_before.tsv
    SRR492065_interleaved.qc.fq.gz
    SRR492065_report.html
    SRR492065_unpaired.qc.fq.gz
    SRR492065_unpaired_summary.tsv

    output/SRR492183/1/qc/0.4.0/fastp:
    SRR492183_fastp.json
    SRR492183_fastp_summary_after.tsv
    SRR492183_fastp_summary_before.tsv
    SRR492183_interleaved.qc.fq.gz
    SRR492183_report.html
    SRR492183_unpaired.qc.fq.gz
    SRR492183_unpaired_summary.tsv
    ```

The fastp output for your data will look like the [following](https://s3.bi.denbi.de/cmg/mgcourses/mg2025/qc/SRR492065_report.html){:target="_blank"}.

!!! question "Task 4"

    How many reads of the wastewater dataset passed the filter?
    Between which quality score ranges does the read1 FASTQ file fall before and after QC.

    ??? Solution
        90.374006% of the reads passed the filter. The quality score for read1 ranges between 26 and 38 before filtering and after QC between 34 and 39. 

#### Host DNA Removal

The Toolkit provides the SRA human-scrubber tool that uses a k-mer based approach to search for human sequences against a k-mer database using human reference sequences.

!!! question "Task 5"

    You can find the result of the tool in the **filterHuman** directory:
    The output of the tool are the filtered sequences, and it also includes statistics about the sequences **before** and **after** filtering.   

    ```BASH
    ls output/*/1/qc/*/filterHuman/
    ```

    The following command tells you how many sequences are left **after** filtering:

    ```BASH
    column -t  -s$'\t' output/*/1/qc/*/filterHuman/*_interleaved_summary_*
    ```

    The **num_seqs** column tells you the number of sequences.

    How many reads were removed from the respective dataset after filtering?

    ??? Solution 
        3768 sequences were detected for SRR492065 and 644 for SRR492183.  

TODO: mention that there are many reads that could be identified ... 

---

### Assembly Results

Once your raw reads have been quality-controlled, we will now perform an assembly using the Metagenomics-Toolkit.
The Toolkit supports several assemblers, one of which is the MEGAHIT assembler, which we will run in this section.

!!! question "Task 6"

    Which tool is used for the assembly? What additional parameters are used?

    ??? Solution 
        MEGAHIT with the additional parameters: minimum contig length of `500 bp` and the preset `meta-sensitive`.
        `meta-sensitive` sets the following parameters: `--min-count 1 --k-list 21,29,39,49,...,129,141`, which causes 
        MEGAHIT to use a longer list of k-mers.
        MEGAHIT is a single-node assembler for large and complex metagenomics
        NGS reads, such as soil. It makes use of a succinct de Bruijn graph
        (SdBG) to achieve a low memory assembly. See the
        [MEGAHIT home page](https://github.com/voutcn/megahit/) for more
        info.


We will now have a first look at some assembly statistics. First of all, locate your assembly results somewhere in your `output` directory.

!!! question "Task 7"
    Where are the assembly results of MEGAHIT stored? And what files are generated per sample?
    
    ??? Solution 
        The assembly results are stored in `output/data/1/assembly/1.2.1/megahit/`
        ```
        cd /vol/volume/sessions/metagenomics_metagenomics-tk
        ```
        ```BASH
        ls -1 output/SRR492065/1/assembly/1.2.3/megahit/
        ```
        ```BASH
            SRR492065_contigs.fa.gz  # the assembled sequences as gzipped fasta
            SRR492065_contigs.fastg  # the assembly graph for inspection for example with Bandage
            SRR492065_contigs_stats.tsv  # some assembly statistics
        ```

!!! question "Task 8"
    Have a look at the `SRR492065_contigs_stats.tsv` file - how large is the assembly and what is the N50?
    
    ??? Solution
        ```BASH
        cd /vol/volume/sessions/metagenomics_metagenomics-tk
        ```
        ```BASH
        cat output/SRR492065/1/assembly/1.2.3/megahit/SRR492065_contigs_stats.tsv | column -s$'\t' -t
        ```
        ```BASH
        SAMPLE     file                     format  type  num_seqs  sum_len   min_len  avg_len  max_len  Q1      Q2      Q3      sum_gap  N50    Q20(%)  Q30(%)  GC(%)
        SRR492065  SRR492065_contigs.fa.gz  FASTA   DNA   1577      11206329  1000     7106.1   245235   1554.0  2659.0  6398.0  0        16801  0.00    0.00    40.46
        ```
        In this assembly run, the assembly length is  11,206,329 bp and the N50 is 16,801 bp. 
 

### Binning Results

Binning is a critical process in metagenomics that involves grouping DNA fragments, known as contigs,
into bins that likely originate from the same organism or genome.
This technique is essential for reconstructing genomes from mixed microbial communities,
where DNA from various organisms is intermingled, and there is often no prior knowledge of the species present.

#### MetaBAT

MetaBAT is an automated metagenome binning software
which integrates empirical probabilistic distances of genome abundance
and tetranucleotide frequency. See the [MetaBAT home page](https://bitbucket.org/berkeleylab/metabat) for more information.

Let's now inspect the bins created by MetaBAT:

!!! question "Task 9"
    How many bins did metabat generate? Locate the metabat results and the fasta files for each bin in the output folder.
    ??? Solution
        ```BASH
        cd /vol/volume/sessions/metagenomics_metagenomics-tk
        ```
        ```BASH
        ls -1 output/SRR492065/1/binning/0.5.0/metabat/
        ```
        ```BASH
        SRR492065_bin.1.fa
        SRR492065_bin.2.fa
        SRR492065_bin.3.fa
        SRR492065_bin.4.fa
        SRR492065_bin.5.fa
        SRR492065_bin_contig_mapping.tsv
        SRR492065_bins_stats.tsv
        SRR492065_contigs_depth.tsv
        SRR492065_notBinned.fa
        ```
        MetaBAT generated 5 bins. The file `SRR492065_notBinned.fa` contains contigs that could not be binned.

The next question you might want to ask is whether you can trust these bins and which organism they represent according to a taxonomy.
Before that, let's have a look at what other information the Metagenomics-Toolkit provides as part of the binning output.

!!! question "Task 10"
    There is a file containing some statistics on the generated bins `SRR492065_bins_stats.tsv`. 
    Find out, which bin has the highest coverage and which one has the highest N50.
    ??? Solution    
        You can have a look on some bin statistics with:
        ```
        cd /vol/volume/sessions/metagenomics_metagenomics-tk
        ```
        ```BASH
        less output/SRR492065/1/binning/0.5.0/metabat/SRR492065_bins_stats.tsv
        ```
        ```BASH
        SAMPLE     BIN_ID              format  type  num_seqs  sum_len  min_len  avg_len  max_len  Q1       Q2       Q3        sum_gap  N50     Q20(%)  Q30(%)  GC(%)  COVERAGE
        SRR492065  SRR492065_bin.1.fa  FASTA   DNA   20        1280335  3167     64016.8  238570   23905.5  45943.5  65927.0   0        185023  0.00    0.00    38.14  169.393
        SRR492065  SRR492065_bin.2.fa  FASTA   DNA   174       2345130  2626     13477.8  68479    6239.0   9839.5   16573.0   0        16802   0.00    0.00    63.56  23.2949
        SRR492065  SRR492065_bin.3.fa  FASTA   DNA   96        1522358  2582     15857.9  63297    6157.0   10569.0  23635.5   0        25145   0.00    0.00    29.95  70.5656
        SRR492065  SRR492065_bin.4.fa  FASTA   DNA   16        1536873  4900     96054.6  245235   22104.0  83532.5  158838.0  0        161674  0.00    0.00    36.70  151.714
        SRR492065  SRR492065_bin.5.fa  FASTA   DNA   436       2720223  2505     6239.0   36553    3386.0   4642.5   7408.5    0        7308    0.00    0.00    32.98  14.7728
        ```
        In this result, bin 1 has the highest coverage (169.393) and highest N50 (185023).
---

➡️ [**Continue to: Assessing Bin Quality**](./bin_quality.md)

