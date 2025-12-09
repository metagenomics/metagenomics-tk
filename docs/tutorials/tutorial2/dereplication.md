In this final part of the tutorial we **aggregate** the MAGs recovered from multiple samples, **dereplicate** them, and briefly inspect both the resulting clusters and the **read-mapping depth** (coverage) of some MAGs in the small example dataset.

Dereplication reduces a large set of MAGs to a **non-redundant genome catalogue** by grouping nearly identical genomes (typically using Average Nucleotide Identity, ANI) and picking a single, 
high-quality representative per cluster. This avoids double-counting strains, simplifies downstream analyses and prevents ambiguous read mapping to many almost-identical genomes.

Within the Metagenomics-Toolkit, dereplication is implemented as a separate module that:

* takes the previous run or a table describing all MAGs (including completeness, contamination and coverage),
* clusters similar genomes in a bottom-up fashion (using Mash and ANI distances), and
* selects the best representative per cluster based on quality and assembly statistics.

---

## Metagenomics-Toolkit Execution

The following lines represent the part of the configuration that tells the Toolkit to run the **Dereplication** module on your existing results.

```YAML linenums="1" title="Dereplication Configuration File Snippet"
---8<--- "default/tutorials/tutorial2/fullPipeline_dereplication.yml:4:41"
```

For reference, your complete parameter file for this final step looks like this:

??? Parameter-file

    ```BASH
    ---8<--- "default/tutorials/tutorial2/fullPipeline_dereplication.yml"
    ```

We assume that all required steps by the previous tutorials are done.

If you run the Toolkit again with the `-resume` flag and this dereplication configuration, it will **reuse** all previous results (QC, assembly, binning, MagAttributes) and only execute the new dereplication step.

!!! Question "Task 1"

    Change into the directory of the Toolkit session.
    ```BASH
    cd /vol/volume/sessions/metagenomics_metagenomics-tk
    ```

    Run the dereplication module on your results using the following command:
    ```BASH
    ---8<--- "scripts/tutorials/tutorial2/dereplication.sh:5:12"
    ```

    This will:

    * collect all the necessary results from your MAGs to create the summary table needed as input,
    * cluster similar genomes based on ANI, and
    * write a `clusters.tsv` file describing the dereplicated genome clusters.

---

## Output

### Dereplication Clusters

The example dereplication input table looks like this (one MAG per row):

```TSV
BIN_ID  COMPLETENESS    COVERAGE        CONTAMINATION   HETEROGENEITY   PATH    N50
SRR492065_bin.3.fa      36.29   171.259 0.39    0       metagenomics-tk/work/58/267ffbad749945cb7dfd98b23a7449/SRR492065_bin.3.fa      73239
SRR492183_bin.2.fa      100.0   142.023 0.1     0       metagenomics-tk/work/90/b5b80ef93e064e7d723eca9bd92241/SRR492183_bin.2.fa      175163
SRR492065_bin.1.fa      85.99   14.7729 18.95   0       metagenomics-tk/work/58/267ffbad749945cb7dfd98b23a7449/SRR492065_bin.1.fa      7308
SRR492065_bin.5.fa      64.75   152.837 0.54    0       metagenomics-tk/work/58/267ffbad749945cb7dfd98b23a7449/SRR492065_bin.5.fa      161674
SRR492065_bin.4.fa      96.32   70.5657 0.14    0       metagenomics-tk/work/58/267ffbad749945cb7dfd98b23a7449/SRR492065_bin.4.fa      25145
SRR492065_bin.2.fa      94.64   23.2949 0.11    0       metagenomics-tk/work/58/267ffbad749945cb7dfd98b23a7449/SRR492065_bin.2.fa      16802
SRR492183_bin.3.fa      92.84   27.7866 0.08    0       metagenomics-tk/work/90/b5b80ef93e064e7d723eca9bd92241/SRR492183_bin.3.fa      21121
SRR492183_bin.1.fa      99.96   45.5264 0.01    0       metagenomics-tk/work/90/b5b80ef93e064e7d723eca9bd92241/SRR492183_bin.1.fa      46795
```

The dereplication module uses these columns as follows:

* `COMPLETENESS`, `CONTAMINATION` – to filter MAGs before clustering.
* `COVERAGE`, `N50`, `HETEROGENEITY` – to decide which genome becomes the **representative** of each cluster.

After running Task 1, dereplication produces a `clusters.tsv` file with one row per MAG and three columns:

* `CLUSTER` – cluster ID (all genomes in the same cluster are near-identical),
* `GENOME` – path or identifier of the genome,
* `REPRESENTATIVE` – `1` if the genome is the chosen representative, `0` otherwise.

A typical file will look like this:

```TSV
CLUSTER  GENOME                                 REPRESENTATIVE
1        bin.1.fa                1
1        bin.2.fa                0
2        bin.8.fasta             1
2        bin.9.fasta             0
3        bin.10.fasta            0
3        bin.32.fa               1
```

*(The exact clustering in your run may differ; this is just an illustrative example.)*

!!! Question "Task 2"

    Locate the dereplication clusters file and inspect how many clusters were formed and which genomes were chosen as representatives.
    Hint: There is a new directory in your output folder.

    ??? Solution
        Change into the Toolkit session directory (if you are not already there):

        ```BASH
        cd /vol/volume/sessions/metagenomics_metagenomics-tk
        ```

        First, locate the new directory in your output folder:

        ```BASH
        ls -la output/
        ```
        There is a new AGGREGATED directory, lets look at our `clusters.tsv` in the similarly named directory.

        ```BASH
        cat output/AGGREGATED/1/dereplication/0.1.1/bottomUpClustering/clusters/clusters.tsv
        ```

        It should look similar to this:

        ```TSV
        CLUSTER GENOME  REPRESENTATIVE
        1       SRR492065_bin.2.fa      0.0
        3       SRR492065_bin.4.fa      1.0
        2       SRR492065_bin.5.fa      0.0
        2       SRR492183_bin.1.fa      1.0
        4       SRR492183_bin.2.fa      1.0
        1       SRR492183_bin.3.fa      1.0
        ```

        This shows you:
        * how many **clusters** exist (number of unique values in the `CLUSTER` column), and  
        * which MAG(s) per cluster have `REPRESENTATIVE=1`. Two MAGs where redundant and are now part of clusters represented only by its best MAG.

---

### Read Mapping Depth (Coverage) of Example MAGs

Coverage values are derived from mapping reads from all input samples back to each representative MAG.

!!! Question "Task 3"

    Let's have a look to examine the coverage (read mapping depth) of our representative MAGs.
    Explore the new directory, find the `readMapping` folder and look for the relative abundance coverage of our representatives with all SRR492065 reads. 

    ??? Solution
        Change again into the Toolkit session directory:

        ```BASH
        cd /vol/volume/sessions/metagenomics_metagenomics-tk
        ```

        Explore the AGGREGATED directory:

        ```BASH
        ls -la output/AGGREGATED/1/
        ```

        There is a `readMapping` directory.
        If we explore further:

        ```BASH
        ls -la output/AGGREGATED/1/readMapping/0.1.0/
        ```
        Also a `genomeCoverage` one.
        If we look at this, we see something like:

        ```
        -rw-r--r-- 1 ubuntu ubuntu  111 Dec  9 18:50 SRR492065_count.tsv
        -rw-r--r-- 1 ubuntu ubuntu  120 Dec  9 18:50 SRR492065_mean.tsv
        -rw-r--r-- 1 ubuntu ubuntu  138 Dec  9 18:50 SRR492065_relative_abundance.tsv
        -rw-r--r-- 1 ubuntu ubuntu  120 Dec  9 18:50 SRR492065_rpkm.tsv
        -rw-r--r-- 1 ubuntu ubuntu  116 Dec  9 18:50 SRR492065_tpm.tsv
        -rw-r--r-- 1 ubuntu ubuntu  118 Dec  9 18:50 SRR492065_trimmed_mean.tsv
        -rw-r--r-- 1 ubuntu ubuntu  110 Dec  9 18:50 SRR492183_count.tsv
        -rw-r--r-- 1 ubuntu ubuntu  113 Dec  9 18:50 SRR492183_mean.tsv
        -rw-r--r-- 1 ubuntu ubuntu  130 Dec  9 18:50 SRR492183_relative_abundance.tsv
        -rw-r--r-- 1 ubuntu ubuntu  110 Dec  9 18:50 SRR492183_rpkm.tsv
        -rw-r--r-- 1 ubuntu ubuntu  112 Dec  9 18:50 SRR492183_tpm.tsv
        -rw-r--r-- 1 ubuntu ubuntu  113 Dec  9 18:50 SRR492183_trimmed_mean.tsv
        drwxr-xr-x 2 ubuntu ubuntu 4096 Dec  9 18:50 logs
        ```
        In `SRR492065_relative_abundance.tsv` we find the right abundances.

        ```BASH
        cat output/AGGREGATED/1/readMapping/0.1.0/genomeCoverage/SRR492065_relative_abundance.tsv
        ```

        ```TSV
        Genome  SRR492065
        unmapped        21.920425
        SRR492065_bin.4 22.608734
        SRR492183_bin.1 45.45554
        SRR492183_bin.2 2.6764407
        SRR492183_bin.3 7.338865
        ```
