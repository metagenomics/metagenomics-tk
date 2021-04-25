

# BWA Workflow

```
./nextflow main.nf -profile slurm -resume --input /raid/meta/test_input.fasta --mapping_samples /raid/meta/test_samples.tsv --list_of_representatives /raid/meta/test/representatives_test.tsv  -entry run_bwa

# Assembly Binning Workflow

```
./nextflow run main.nf -profile slurm -resume  --interleaved --input /raid/meta/test/assembly_binning/test_samples_interleaved.tsv --megahit --metabat --getreads --checkm --database /vol/spool/CAMI/CAMI_MOUSEGUT/reassembly/checkmdb  --ending .fa -entry run_assembly_binning
```
