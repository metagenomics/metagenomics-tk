# Run Fragment Recruitment

The fragment recruitment module can be used to find genomes in a set of read datasets.

In case the fragment recruitment module is part of the full pipeline pr per-sample pipeline configuration then
reads that could not be mapped back to a contig are screened for a user provided list of MAGs.
Detected genomes are included in all other parts of the remaining pipeline.
Look out for their specific headers to differentiate results based on real assembled genomes and the reference genomes.

**Note:** This module currently only supports illumina data. 

## Input

=== "Command"

    ```
    -entry wFragmentRecruitment -params-file example_params/fragmentRecruitment.yml
    ```

=== "Configuration file for fragment recruitment via mash screen and BWA"

    !!! warning "Warning"
     
        **The configuration file shown here is for demonstration and testing purposes only. 
          Parameters that should be used in production can be viewed in the fragment recruitment section 
          of one of the yaml files located in the `default` folder of the Toolkit's Github repository.**

    ```YAML
    ---8<--- "example_params/fragmentRecruitmentMashScreen.yml"
    ``` 

=== "Configuration file for fragment recruitment via mash screen and Bowtie"


    !!! warning "Warning"
     
        **The configuration file shown here is for demonstration and testing purposes only. 
          Parameters that should be used in production can be viewed in the fragment recruitment section 
          of one of the yaml files located in the `default` folder of the Toolkit's Github repository.**

    ```YAML
    ---8<--- "example_params/fragmentRecruitment_fraction/fragmentRecruitmentMashScreenBowtie.yml"
    ``` 

=== "Input TSV file for genomes"

    ```TSV
    ---8<--- "test_data/fragmentRecruitment/mags.tsv"
    ```

=== "Input TSV file for paired end reads"

    ```TSV
    ---8<--- "test_data/fragmentRecruitment/paired.tsv"
    ```

=== "Input TSV file for single end reads"

    ```TSV
    ---8<--- "test_data/fragmentRecruitment/single.tsv"
    ```

**NOTE!** The file names of all provided genomes must be unique.

The following parameters can be configured:

  * mashDistCutoff: All hits below this threshold are  discarded.

  * mashHashCutoff: All hits that have a lower count of matched minimum hashes are discarded.

  * coveredBasesCutoff: Number of bases that must be covered by at least one read. By how many reads
    the bases must be covered can be configured via the coverm setting (coverm: "  --min-covered-fraction 0  ").

## Output

The module outputs mash screen and bowtie alignment statistics. 
Furthermore, the module provides a coverm output which basically reports all metrics
about the found genomes (e.g covered bases,length, tpm, ...).
