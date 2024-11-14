This section describes initial steps to run the pipeline on a cloud-based cluster system.
Other parameters are found in the configuration part.

## Full Pipeline vs. Seperate Modules

The Metagenomics-Toolkit allows you to run either the full pipeline of assembly, binning and many other downstream analysis tasks or individual modules.
The toolkit can be configured by providing the module configuration via a yml file and a flag for the corresponding [module](modules/introduction.md) or [full pipeline mode](full_pipeline.md). Options for the global pipeline configuration can be viewed [here](configuration.md).

The Full Pipeline mode consists mainly of two parts. One part of the Toolkit processes each dataset individually by applying modules such as quality control, assembly, binning and annotation. The second part aggregates and combines the results of the individual datasets by applying modules such as dereplication and
co-occurrence analysis. You can either run the first and second part consecutively or seperately where the second part can be applyed on the output of first 
part.

## Error Strategy 

All tools follow the same error strategy. The execution of a tool is retried three times. If the run fails the fourth time, it will be ignored.
If the execution is ignored, the toolkit will continue to run all tools that do not depend on the output of the failed tool run.
Exceptions of this handling are specified in the corresponding module section.

*Note!* Please do never place sensitive information in any of the yml configuration files since the configuration is part of the pipeline output.

