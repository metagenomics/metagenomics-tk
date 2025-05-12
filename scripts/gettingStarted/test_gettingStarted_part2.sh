# Download files to demonstrating the local path case
mkdir inputFiles
wget -P inputFiles https://openstack.cebitec.uni-bielefeld.de:8080/swift/v1/meta_test/small/read1_1.fq.gz 
wget -P inputFiles https://openstack.cebitec.uni-bielefeld.de:8080/swift/v1/meta_test/small/read2_1.fq.gz 

# Create a file that describes your input data 
INPUT_FILES=$(pwd)/inputFiles/input.tsv
echo -e "SAMPLE\tREADS1\tREADS2" > $INPUT_FILES
echo -e "MYDATA\t$(readlink -f inputFiles/read1*)\t$(readlink -f inputFiles/read2*)" >> $INPUT_FILES

# Run the following command with files stored on local disk
NXF_HOME=$PWD/.nextflow NXF_VER=25.04.0 nextflow run metagenomics/metagenomics-tk -work-dir $(pwd)/work \
    -profile slurm \
    -entry wFullPipeline \
    -ansi-log false \
    -params-file  https://raw.githubusercontent.com/metagenomics/metagenomics-tk/refs/heads/master/default/fullPipeline_illumina_nanpore_getting_started_part1.yml \
    --logDir log2 \
    --s3SignIn false \
    --scratch /vol/scratch \
    --databases /vol/scratch/databases \
    --output my_data_output \
    --input.paired.path inputFiles/input.tsv


make check LOG_DIR=$(pwd)/log2
