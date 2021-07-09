## Output and best practice

### Motivation

* The output section is a collection of `best practices` for storing results of the `meta-omics-toolkit` output.
The definitions are motivated by the fact that the pipeline will be continuously updated and results of different pipeline
versions and modes must be differentiated.

* The idea is to run the pipeline on results of previous runs.

### Rules for dataset output

Outputs are produced by using the `publish dir` directive.

```
DATASET_ID/RUN_ID/MODULE/VERSION/TOOL/
```
where
   * `DATASET_ID` specifies the ID of a dataset such as the SRA dataset.
   * `RUN_ID` specifies one possible run of the full or partial pipeline. The `RUN_ID` identifier can be any user provided identifier to keep track of multiple pipeline runs.
   * `MODULE` specifies the name of the pipeline module (e.g. binning).
   * `VERSION` specifies the module version number which follows semantic versioning (1.2.0).
   * `TOOL` specifies the name of the tool that is executed as part of the module (e.g `megahit` of the assembly module).

It is suggested that a RUN_ID output should never contain multiple versions of the same module. E.g.: 
`DATASET_ID/1/Binning/1.2.0/metabat` and `DATASET_ID/1/Binning/1.3.0/metabat`.

If a partial pipeline run (B) uses outputs of a previous run (A) (e.g. a binning tool uses the output of an assembler) and the previous run (A) alreads contains
the output of an older version of run (B), then a new RUN_ID folder must be created.

If a partial pipeline run (B) uses outputs of a previous run (A) (e.g. a binning tool uses the output of an assembler) and the previous run (A) **does not** contain
the output of an older version of run (B), then the existing RUN_ID folder must be **reused**.

#### Run Versioning

Every dataset must contain a `TOOL` folder called `config`. The `config` folder contains descriptions of the parameters and the version used for the specific pipeline run.

### Examples

## Example 1:

We assume that the following folder already exists:

```
/SRA1/1/ASSEMBLY/1.2/MEGAHIT
```


If the MODULE output does not contain a BINNING output then the existing RUN folder must be reused:

```
/SRA1/1/ASSEMBLY/1.2/MEGAHIT
/SRA1/1/BINNING/0.3/METABAT
```

## Example 2:

We assume that the following folders already exists:

```
/SRA1/1/ASSEMBLY/1.2/MEGAHIT
/SRA1/1/BINNING/0.3/METABAT
```

If the MODULE output does contain a BINNING output then a new RUN folder must be created:

```
/SRA1/1/ASSEMBLY/1.2/MEGAHIT
/SRA1/1/BINNING/0.3/METABAT
/SRA1/2/BINNING/0.4/METABAT
```
