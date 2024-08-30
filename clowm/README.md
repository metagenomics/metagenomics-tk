# Meta-Omics-Toolkit

The Meta-Omics-Toolkit allows you to run either the full pipeline of assembly, binning and many other downstream analysis
tasks or individual modules. The toolkit can be configured by providing the module configuration via a yml file and a flag
for the corresponding module or full pipeline mode.

All tools follow the same error strategy. The execution of a tool is retried three times. If the run fails the fourth time,
it will be ignored. If the execution is ignored, the toolkit will continue to run all tools that do not depend on the output
of the failed tool run. Exceptions of this handling are specified in the corresponding module section.
