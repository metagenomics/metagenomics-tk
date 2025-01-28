## Contents

The Getting Started tutorial explains the first steps in running the Toolkit on a cloud-based cluster system.
It consists of two main parts. 

1. In the [first part](full_pipeline.md) you will learn how to run and configure the Toolkit.
2. The [second part](aggregation.md) tells you how to aggregate multiple samples of a Toolkit output.

Other parameters can be found in the [configuration](configuration.md) section.

## Full Pipeline vs. Seperate Modules

The Metagenomics-Toolkit allows you to run either the full pipeline of assembly, binning and many other downstream analysis tasks or the individual modules.
The Toolkit can be configured by providing the module configuration via a yml file and a flag for the corresponding [module](modules/introduction.md) or [full pipeline mode](full_pipeline.md). Options for the global pipeline configuration can be viewed [here](configuration.md).

The Full Pipeline mode consists mainly of two parts. One part (per-sample) of the Toolkit processes each dataset individually 
by applying analyses such as quality control, assembly, binning and annotation. The second part aggregates and combines the results of the individual datasets by applying modules such as dereplication and
co-occurrence analysis. You can either run the first and second part consecutively or seperately where the second part can be applied on the output of first 
part.

!!! note "Sensitive data"
    Please do never place sensitive information in any of the yml configuration files since the configuration is part of the pipeline output.

### Error Strategy 

All tools follow the same error strategy. The execution of a tool is retried three times. If the run fails the fourth time, it will be ignored.
If the execution is ignored, the Toolkit will continue to run all tools that do not depend on the output of the failed tool run.
Exceptions of this handling are specified in the corresponding module section.
