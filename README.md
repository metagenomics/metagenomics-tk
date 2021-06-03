
# Pipeline Modes

## Run Full Pipeline

```
./nextflow run main.nf -work-dir /shared/directory/test -profile PROFILE  -resume -entry run_ppeline -params-file example_params/full_pipeline_params.yml
```

where
 *  /shared/directory/test is a directory that is shared between multiple machines.
 * PROFILE can be either `local,conda` or `slurm,conda` depending on which environment the pipeline should be executed.

## BWA Workflow

```
./nextflow main.nf -profile slurm -resume --input /raid/meta/test_input.fasta --mapping_samples /raid/meta/test_samples.tsv --list_of_representatives /raid/meta/test/representatives_test.tsv  -entry run_bwa
```

## Assembly Binning Workflow

```
./nextflow run main.nf -profile slurm -resume  --interleaved --input /raid/meta/test/assembly_binning/test_samples_interleaved.tsv --megahit --metabat  --checkm --database /vol/spool/CAMI/CAMI_MOUSEGUT/reassembly/checkmdb  --gtdb_database /vol/spool/CAMI/CAMI_MOUSEGUT/methods/segata/bowtie/database/release89 --gtdb --ending .fa --getreads  --output=/raid/meta/out  -entry run_assembly_binning -with-timeline timeline.html -with-trace -with-report report.html
```

## Metabolic Reconstruction

An environment must be created using 

```
conda env create -n metabolic -f conda/metabolic.yml
```

You can also update an existing environment using:

```
conda env update --n metabolic --file conda/metabolic.yml
```

Further you have to install the python packaes

```
pip install -r requirements/metabolic.txt
```

Install IBM cplex 

You can run the pipeline by specifying the environment:

```
./nextflow  /vol/spool/meta/main.nf -profile conda -resume  --environment /vol/spool/conda/anaconda3/envs/metabolic  --input /vol/spool/meta/test/data/table.tsv
```

# Developer Guidelines

## Workflows

1. Worfklow names that can not be used directly and are just meant for internal use should start with an underscore.

2. Every Workflow should contain a short description of the functionality. 

3. Workflow names must start with `w`. 

## Process

1. Process names should start `p`

## Other

1. Magic numbers hould not be used.

2. Variable, method, workflow, folder and process names should be written in camelcase.
