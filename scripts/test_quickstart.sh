# This file will be referenced on the online wiki

git clone git@github.com:metagenomics/metagenomics-tk.git
cd metagenomics-tk

make nextflow

./nextflow run main.nf \
	  -profile standard \
	  -entry wFullPipeline \
	  -params-file default/quickstart.yml \
	  --s3SignIn false \
	  --scratch false \
	  --output output \
	  --databases $(pwd)/databases
