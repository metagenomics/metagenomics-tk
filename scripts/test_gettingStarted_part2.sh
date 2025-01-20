# Download files to demonstrating the local path case
mkdir inputFiles
wget -P inputFiles https://openstack.cebitec.uni-bielefeld.de:8080/swift/v1/meta_test/small/read1_1.fq.gz 
wget -P inputFiles https://openstack.cebitec.uni-bielefeld.de:8080/swift/v1/meta_test/small/read2_1.fq.gz 

# Create a file that describes your input data 
INPUT_FILE=$(pwd)/inputFiles/input.tsv
echo "MYDATA\tREADS1\tREADS2" > $INPUT_FILES
echo "MYDATA\t$(readlink -f inputFiles/read1*)\t$(readlink -f inputFiles/read2*)" >> $INPUT_FILES

# Run the following command with files stored on local disk
NXF_HOME=$PWD/.nextflow NXF_VER=23.10.0 nextflow run main.nf -work-dir $(pwd)/work \
    -profile slurm \
    -entry wFullPipeline \
    -params-file https://github.com/metagenomics/metagenomics-tk/blob/master/default/fullPipeline_illumina_nanpore_getting_started_part1.yml
    --s3SignIn false \
    --scratch /vol/scratch \
    --databases /vol/scratch/databases \
    --output my_data_output \
    --input.paired.path inputFiles/input.tsv

