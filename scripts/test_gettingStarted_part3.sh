# Run the following command to test the MetaSpades assembler 
NXF_HOME=$PWD/.nextflow NXF_VER=23.10.0 nextflow run metagenomics/metagenomics-tk -work-dir $(pwd)/work \
    -profile slurm \
    -entry wFullPipeline \
    -params-file  https://raw.githubusercontent.com/metagenomics/metagenomics-tk/refs/heads/feat/getting_started_guide/default/fullPipeline_illumina_nanpore_getting_started_part3.yml \
    --s3SignIn false \
    --scratch /vol/scratch \
    --databases /vol/scratch/databases \
    --output my_data_spades_output \
    --input.paired.path inputFiles/input.tsv

cat my_data_spades_output/MYDATA/1/assembly/*/metaspades/MYDATA_contigs_stats.tsv
SAMPLE  file                    format  type    num_seqs  sum_len    min_len avg_len max_len Q1      Q2      Q3      sum_gap N50     Q20(%)  Q30(%)  GC(%)
MYDATA  MYDATA_contigs.fa.gz    FASTA   DNA     1761      4521427    1000    2567.5  26470   1149.0  1408.0  2119.0  0       3799    0.00    0.00    55.07

cat my_data_output/MYDATA/1/assembly/*/megahit/MYDATA_contigs_stats.tsv
SAMPLE  file                    format  type    num_seqs  sum_len    min_len avg_len max_len Q1      Q2      Q3      sum_gap N50     Q20(%)  Q30(%)  GC(%)
MYDATA  MYDATA_contigs.fa.gz    FASTA   DNA     95227     60110517   56      631.2   346664  234.0   286.0   439.0   0       1229    0.00    0.00    57.34
