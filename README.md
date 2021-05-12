

# BWA Workflow

```
./nextflow main.nf -profile slurm -resume --input /raid/meta/test_input.fasta --mapping_samples /raid/meta/test_samples.tsv --list_of_representatives /raid/meta/test/representatives_test.tsv  -entry run_bwa
```

# Assembly Binning Workflow

```
./nextflow run main.nf -profile slurm -resume  --interleaved --input /raid/meta/test/assembly_binning/test_samples_interleaved.tsv --megahit --metabat  --checkm --database /vol/spool/CAMI/CAMI_MOUSEGUT/reassembly/checkmdb  --gtdb_database /vol/spool/CAMI/CAMI_MOUSEGUT/methods/segata/bowtie/database/release89 --gtdb --ending .fa --getreads  --output=/raid/meta/out  -entry run_assembly_binning -with-timeline timeline.html -with-trace -with-report report.html
```

# Metabolic Reconstruction

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

# Run Full Pipeline

```
./nextflow run main.nf -work-dir /raid1/test -profile local,conda -resume -entry run_pipeline -params-file full_pipeline_params.yml
```
