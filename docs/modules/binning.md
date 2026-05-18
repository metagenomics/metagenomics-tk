# Binning 

The binning module groups together assembled contigs to create bins representing genomes.
Depending on whether ONT, short read and/or multi-sample binning should be executed, 
multiple binning configurations are offered.

For multi-sample binning, a group of samples must first be configured. 
If one of these samples fails QC or assembly, leaving only one sample per group,
then the sample is switched from multi-sample binning to per-sample binning.

## Short Read 

For short read data the Toolkit supports MetaBAT2 and Semibin2.

### Input

=== "Command for short read data with optional single end reads (Metabat2)"

    ```
    ---8<--- "scripts/modules/binning/test_binning.sh:func"
    ```

=== "TSV table short read"

    ```TSV
    ---8<--- "test_data/binning/samples.tsv"
    ```

=== "TSV table unpaired short read"

    ```TSV
    ---8<--- "test_data/binning/samplesUnpaired.tsv"
    ```

=== "TSV table for contigs"

    ```TSV
    ---8<--- "test_data/binning/assembly.tsv"
    ```

### Output

#### Contig Coverage

CoverM is executed on all alginment files to provide basic contig coverage information.
Here the output are two files. One file is the output of all coverm default methods.
The other file are metabat specific metrics which are offered by CoverM.

`contigCoverage/*_default_coverm_coverage.tsv`

: The following methods are used as parameters for CoverM to produce this file `mean trimmed_mean variance length count reads_per_base rpkm tpm`.

`contigCoverage/*_metabat_coverm_coverage.tsv`

: The parameter for CoverM is just `metabat`. 

#### Alignment

`contigMapping/*.bam`

: Alignment file containing the alignment of reads mapped back to the assembly.

`contigMapping/*_unmapped.fq.gz`

: Reads that could not be mapped back.

#### Genome Coverage

CoverM is executed on all alginment files to provide basic coverage information of all generated genomes.
Every file is based on a CoverM method `mean trimmed_mean variance length count reads_per_base rpkm tpm`.

`genomeCoverage/*_count.tsv`

`genomeCoverage/*_mean.tsv`

`genomeCoverage/*_relative_abundance.tsv`

`genomeCoverage/*_rpkm.tsv`

`genomeCoverage/*_tpm.tsv`

`genomeCoverage/*_trimmed_mean.tsv`

#### Read Mapping Quality

`readMappingQuality/*_flagstat.tsv`

: Samtools flagstat output

`readMappingQuality/*_flagstat_failed.tsv`

: Samtools flagstat output of failed reads in horizontal format.

`readMappingQuality/*_flagstat_passed.tsv`

: Samtools flagstat output of passed reads in horizontal format.

#### Binning Tool Output

The output has always the same directory structure for all binning tools:

Example Semibin2:

`semibin2/*_bin.*.fa`

: Generated bins.

`semibin2/*_bin_contig_mapping.tsv`

: Bin to contig mapping.

`semibin2/*_bins_stats.tsv`

: Basic sequence stats produced by `seqkit`. 

`semibin2/*_contigs_depth.tsv`

: Output of the jgi_summarize_bam_contig_depths file.

`semibin2/*_notBinned.fa`

: Contigs that could not be binned.

## ONT

For Nanopore data the Toolkit supports MetaBAT2 and MetaCoAG.

### Input

=== "Command for Nanopore reads (Metabat2)"

    ```BASH
    ---8<--- "scripts/modules/binning/test_binningONT.sh:func"
    ```

=== "TSV table reads"

    ```TSV
    ---8<--- "test_data/binningONT/samplesONT.tsv"
    ```

=== "TSV table for contigs"

    ```TSV
    ---8<--- "test_data/binningONT/assemblyONT.tsv"
    ```

=== "TSV for quality of ONT reads"

    ```TSV
    ---8<--- "test_data/binningONT/quality.tsv"
    ```

### Output

The output is the same as for short read data. 

## Multi Sample Binning Short Read

For short read multi sample data the Toolkit supports SemiBin2.

### Input

=== "Command for Nanopore reads (SemiBin2)"

    ```BASH
    ---8<--- "scripts/modules/binning/test_multiBinning.sh:func"
    ```

=== "TSV for quality of paired reads"

    ```TSV
    ---8<--- "test_data/multiBinning/samples.tsv"
    ```

=== "TSV for quality of unpaired reads"

    ```TSV
    ---8<--- "test_data/multiBinning/samplesUnpaired.tsv"
    ```

=== "TSV for contigs"

    ```TSV
    ---8<--- "test_data/multiBinning/assembly.tsv"
    ```

=== "TSV for group information"

    ```TSV
    ---8<--- "test_data/multiBinning/groups.tsv"
    ```

### Output

In addition to the files produced for short read data, there are also
the following files created:

#### Concatenated Assembly Alignment

`concatenatedAssemblyMapping/*.bam`

: This is an alignment file containing the alignment of reads from every sample mapped back to the concatenated assembly consisting of all the assemblies specified in the groups file.

`concatenatedAssemblyMapping/*_unmapped.fq.gz`

: Reads that could not be mapped back.

## Multi Sample Binning ONT

### Input

=== "Command for Nanopore reads (SemiBin2)"

    ```BASH
    ---8<--- "scripts/modules/binning/test_multiBinningONT.sh:func"
    ```

=== "TSV for ONT reads"

    ```BASH
    ---8<--- "test_data/multiBinningONT/samplesONT.tsv"
    ```

=== "TSV for contigs"

    ```BASH
    ---8<--- "test_data/multiBinningONT/assemblyONT.tsv"
    ```

=== "TSV for group information"

    ```TSV
    ---8<--- "test_data/multiBinningONT/groups.tsv"
    ```

=== "TSV for quality information"

    ```TSV
    ---8<--- "test_data/multiBinningONT/quality.tsv"
    ```

### Output

The output is the same as for multi sample short read data. 
