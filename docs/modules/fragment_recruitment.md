# Run Fragment Recruitment

The fragment recruitment module can be used to find genomes in a set of read datasets.

## Input

=== "Command"

    ```
    -entry wFragmentRecruitment -params-file example_params/fragmentRecruitment.yml
    ```

=== "Configuration file for fragment recruitment via mash screen"

    ```YAML
    ---8<--- "../../example_params/fragmentRecruitmentMashScreen.yml"
    ``` 

=== "Configuration file for fragment recruitment via frhit (deprecated)"

    ```YAML
    ---8<--- "../../example_params/fragmentRecruitment.yml"
    ```

=== "Input TSV file for genomes"

    ```TSV
    ---8<--- "../../test_data/fragmentRecruitment/mags.tsv"
    ```

=== "Input TSV file for paired end reads"

    ```TSV
    ---8<--- "../../test_data/fragmentRecruitment/paired.tsv"
    ```

=== "Input TSV file for single end reads"

    ```TSV
    ---8<--- "../../test_data/fragmentRecruitment/single.tsv"
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
