# Meta-Omics-Toolkit

## Modules

All module configurations are the same as the full pipeline run with the sole difference that entry and param-file parameters are different.

*Note!* Please do never place sensitive information in any params file since the content is upload as part of the pipeline run.

### Run Full Pipeline

```
./nextflow run main.nf -work-dir /shared/directory/test -profile PROFILE  -resume -entry wPipeline -params-file example_params/fullPipeline.yml
```

where
 *  /shared/directory/test is a directory that is shared between multiple machines.
 * PROFILE can be either `standard` (local use) or `slurm` depending on which environment the pipeline should be executed.

**Note!** Metabolomics part is currently excluded from full pipeline run.

#### Input

* [params-file](example_params/fullPipeline.yml)

* [Tsv Table](test_data/fullPipeline/reads_split.tsv): Must include the columns `SAMPLE`, `READS1` and `READS2`. `SAMPLE` must contain unique dataset identifiers
without whitespaces or special characters. `READS1` and `READS2` are paired reads and can be HTTPS URLs, S3 links or files.

#### Output (Overview)

In addition to the pipeline module outputs defined in the next sections (Dereplication, MagAttributes, etc), the following outputs are produced. 

 * quality control (fastq) 

 * assembly (contigs)

 * binning (genomes)
 
 * read mapping (bam files)

#### Configuration

Options for global pipeline configuration can be viewed [here](docs/configuration/global.md).

##### Additional Configuration

Nextflow usually stores downloaded files in the work directory. If enough scratch space is available on the worker nodes then this can be prevented by specifying
s3 links ([example](test_data/fullPipeline/reads_split_s3.tsv)) in the input tsv file and `download` parameter in the input yaml ([example](example_params/fullPipelineQC.yml)).

### Run Fragment Recruitment

```
-entry wFragmentRecruitment -params-file example_params/fragmentRecruitment.yml
```

### Assembly

```
-entry wAssembly -params-file example_params/assembly.yml
```

#### Input

* [params-file](example_params/assembly.yml)

* [Sample TSV Table](test_data/assembly/samples.tsv)

#### Output

The output is a gzipped fasta file containing contigs.

#### Error Handling

On error with exit codes ([-9, 137]) (e.g. due to memory restrictions), the tool is executed again with higher cpu and memory values.
The memory and cpu values are computed by the formula 2^(number of attempts) * (cpu/memory value of the assigned flavour).
The highest possible cpu/memory value is restricted by the highest cpu/memory value of all flavours defined in the resource section 
(see global [configuration](docs/configuration/global.md) section). 

### Cooccurrence

```
-entry wCooccurrence -params-file example_params/coocurrence.yml
```

#### Input

* [params-file](example_params/coocurrence.yml)

* [abundance TSV table](test_data/assembly/samples.tsv): Contains abundance values of mags per sample.

* [gtdb TSV table](test_data/assembly/samples.tsv): GTDB assignmend of all samples that were produced by magAttributes module.

#### Output

 * Graphml file for further processing 

### Dereplication

```
-entry wDereplication -params-file example_params/dereplication.yml
```

#### Input

* [params-file](example_params/dereplication.yml)

* [Tsv Table](test_data/dereplication/input.tsv): Must include the columns `DATASET`, `BIN_ID`, `PATH`, `COMPLETENESS`, `CONTAMINATION`, `COVERAGE`, `N50` and `HETEROGENEITY`. 
Completeness and contamination can be used for filtering (see `params-file`). `N50`, `COVERAGE` and `HETEROGENEITY` are used for selecting the representative of every cluster.
You can set values of these columns to zero if data is not available or if you don't want the representative selection to be influenced by theses columns. Make sure that `BIN_ID` is a unique identifier.

#### Output

The output tsv file (`final_clusters.tsv`) contains the columns `CLUSTER`, `GENOME` and `REPRESENTATIVE` where `CLUSTER` identifies a group of genomes, `GENOME` represents the path or
link of a genome and `REPRESENTATIVE` is either 0 or 1 (selected as representative).
If `sans` is specified in the configuration file (see examples folder), then [SANS](https://gitlab.ub.uni-bielefeld.de/gi/sans) is used to dereplicate every cluster reported by the previous step further down 
to generate strain specific clusters. 

### Read Mapping

```
-entry wReadMapping -params-file example_params/readMapping.yml
```

#### Input

* [params-file](example_params/readMapping.yml)

* [MAGs TSV Table](test_data/readMapping/mags.tsv): 

* [SAMPLES TSV Table](test_data/readMapping/samples.tsv): 

#### Output

The produced output files are the following: count.tsv, mean.tsv, mean_mincov10.tsv, rpkm.tsv, tpm.tsv, trimmed_mean.tsv.
The content of the files are produced by coverm. All metrics are explained on the coverm github page: https://github.com/wwood/CoverM .

### MagAttributes

```
 -entry wMagAttributes -params-file example_params/magAttributes.yml 
```

#### Input

* [params-file](example_params/magAttributes.yml)

* [Tsv Table](test_data/magAttributes/input.tsv): Must include at least `DATASET` identifier and mag specific `PATH` and `BIN_ID` column.

#### Output


##### Prokka

Prokka computes `*.err`, `*.faa`, `*.ffn`, `*.fna`, `*.fsa`, `*.gbk`, `*.gff`, `*.sqn`, `*.tbl`, `*.tbl` for every bin.
Details of all files can be read on the Prokka page.
In addition it also computes a summary tsv file which adheres to the magAttributes specification.

##### GTDBTk

All GTDB files include the GTDB specific columns in addition to a `SAMPLE` column (`SAMPLE_gtdbtk.bac120.summary.tsv`, `SAMPLE_gtdbtk.ar122.summary.tsv`).
In addition, this module produces a file `SAMPLE_gtdbtk_CHUNK.tsv` that combines both files and adds a `BIN_ID` column that adheres to the magAttributes specification

##### Checkm

The Checkm output adheres to the magAttributes specification and adds a `BIN_ID` and `SAMPLE` column to the output file.

### Plasmids

The plasmid module is able to identify contigs as plasmids and it is able to assemble plasmids from the samples fastq data.
For running the plasmid assembly we suggest to run the full pipeline mode with the enabled [plasmid module](example_params/fullPipeline_fraction/fullPipeline_fraction_plasmid.yml).

```
-entry wPlasmidsPath -params-file example_params/plasmids.yml
```

#### Input

* [params-file](example_params/plasmid.yml)

* [Tsv Table](test_data/plasmid/input_contigs.tsv)

#### Output

##### SCAPP

SCAPP detects plasmid sequences out of the samples assembly graph.
It reports sequences as gzipped fasta files (`*_plasmids.fasta.gz`). A basic statistic (`*_plasmids_stats.tsv`)
is also generated.

##### PlasClass

PlasClass is able to identify contigs as plasmids. It reports gzipped fata files (`*_plasmids.fasta.gz`)
and their probabilities (`*_probs.tsv`).

##### PLSDB

PLSDB includes a curated set of plasmid sequences that were extracted from databases like refseq.
The metadata of found sequences are reported in `*.tsv` and the metadata of the filtered sequences in `*_kmerThreshold_X.tsv`.

### Annotation

```
-entry wAnnotateLocal -params-file example_params/annotation.yml
```

#### Input  

* [params-file](example_params/annotation.yml)

* [Tsv Table](test_data/annotation/input_small.tsv)

#### Output

##### Diamond

Calculated significant matches of a nucleotide/protein query which was compared against a database.

##### Prodigal

Predicted genes in the `*.faa` `*.fna` and `*.gff` file format.

##### KEGGFromDiamond

Result `*.tsv` file filled with KEGG informations (linke modules, KO's, ...) which could be linked to the input Diamond hits.
  
##### Resistance Gene Identifier (rgi)

The `*rgi.tsv` files contain the found CARD genes.

## Error Strategy

All tools follow the same error strategy. The execution of a tool is retried four times. If the run fails the fourth time, it will be ignored.
If the execution is ignored, the toolkit will continue to run all tools that do not depend on the output of the failed tool run.
Exceptions of this handling are specified in the corresponding module section.

## S3 Configuration

All modules of the pipeline can be used in conjunction with S3.
You will have to create a configuration file that can be provided to nextflow with " -c " Parameter.

```
aws {
  accessKey = 'xxx'
  secretKey = 'xxx'

    client {
      s_3_path_style_access = true
      connectionTimeout = 120000
      maxParallelTransfers = 28 
      maxErrorRetry = 10
      protocol = 'HTTPS'
      endpoint = 'https://openstack.cebitec.uni-bielefeld.de:8080'
      signerOverride = 'AWSS3V4SignerType'
    }
}
```

If you want to upload tool results to s3, just update the output parameter in the configuration file from `/path/to/output` to `s3://bucket/path/to/output`

If you want to use the annotation module, you also have to provide your credentials in an aws credential style.
Create a file that looks like this and fill in your credentials:

```
[default]
aws_access_key_id=ABCDEKEY
aws_secret_access_key=ABCDEKEY
```
You have to reference this file in the annotation parameter yml's s5cmd keyfiles section.  

## Other 

* [Developer Guidelines](docs/developer_guidelines.md)

* [Module Specification](docs/specification.md)
