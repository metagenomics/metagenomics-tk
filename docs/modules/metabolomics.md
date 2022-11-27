# Metabolomics

The metabolomics module runs genome scale metabolic modeling analysis based on a supplied genome or directly on its proteins.
The module is able to use gapseq and carveme for analysing genomes and carveme for analysing predicted proteins
which depends on the configuration you provide as input.

**Note:** If carvem is specificed in fullPipeline mode then carveme is executed with proteins as input.

All generated models are used for further downstream analysis such as the "Minimum Resource Overlap" computation by smetana.

## Input

=== "Command"

    ```
    -entry wMetabolomics -params-file example_params/metabolomics
    ```

=== "Configuration file for providing genomes"

    ```YAML
    ---8<--- "../example_params/metabolomics.yml"
    ```

Almost all tools of this module are using linear programming solvers. The tool developers are recommending the use of the cplex solver
that is included in the [IBM ILOG CPLEX Optimization Studio](https://www.ibm.com/de-de/products/ilog-cplex-optimization-studio) which is free for students and academics
through the IBM Academic Initiative programm. 
Since the toolkit uses docker images that are downloaded from public Docker Hub repositories and the cplex license is not allowed
to be distributed, we prepared a [Dockerfile](../cplex/docker/Dockerfile) that allows you to build your own local docker image with all metabolomics specific tools installed.
Just copy your cplex binary to the cplex docker folder and build your own docker image. You can override all existing images via the command line.

In the following example your the image name is metabolomics:0.1.0:

* `--gapseq_image=metabolomics:0.1.0` (Optional)
* `--smetana_image=metabolomics:0.1.0` (Required)
* `--carveme_image=metabolomics:0.1.0` (Required)
* `--memote_image=metabolomics:0.1.0` (Optional. Memote is not able to detect the solver automatically. Please specify `--solver` in the configuration file if you are not using the glpk solver.)

For gapseq and memote we are using a publicly available docker image that uses the freely available glkp solver which means that you don't have to provide this parameter.
If you want to build your own image, please use the `beforeProcessScript` parameter. This parameter expects a bash script that accepts the docker image name as a parameter.
The script is executed right before the actual docker image is started. 
You could for example provide a script that builds the actual image right before running the tool on the VM. 
It would be also possible to push the docker image to a private dockerhub repository and login to your docker account via this script.
We provide two example template scripts in the cplex folder.
Please note that in both cases you distribute the docker image with your cplex binary on all machines where you run the toolkit.
If you login to dockerhub then your credentials will saved on the VM. If you are not the only docker user on the machine we do not recommend this approach! 

## Output

### GapSeq / CarveMe 

Both tools are generating genome scale metabolic reconstruction models (`*.xml`). 
All models are translated to json format and substrats, products and reactions are saved in distinct files.

### Memote

Memote tests metabolic reconstruction models and therefore produces a machine readable json file  (`*_report.json.gz`)
and a human readable tsv (`*_metrics.tsv`) and html (`*_report.html`) file.

### Smetana

Smetana is used for analysing possible interactions in microbial communities. Smetana`s global and detailed modes 
are executed per sample. The Smetana output is saved in `*_detailed.tsv` and `*_global.tsv`.
