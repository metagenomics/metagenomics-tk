We are going to run an assembly with the metagenomics toolkit using MEGAHIT as an assembler, the toolkit is also capable of running metaSPades. 

## MEGAHIT 

MEGAHIT is a single node assembler for large and complex metagenomics
NGS reads, such as soil. It makes use of succinct de Bruijn graph
(SdBG) to achieve low memory assembly. MEGAHIT can optionally utilize
a CUDA-enabled GPU to accelerate its SdBG contstruction. See the
`MEGAHIT home page <https://github.com/voutcn/megahit/>`_ for more
info.


## Metagenomics-Toolkit  


The following lines need to be added to your parameter file in order to run the assembly (below the qc-part):
```BASH
---8<--- "default/tutorials/tutorial1/fullpipeline_assembly.yml:28:38"
```


The following command runs the quality control and assembly module of the Toolkit on one machine.
```BASH
---8<--- "scripts/tutorials/tutorial1/test_assembly.sh:3:12"
```



For reference, your complete parameter file should look like this:
```BASH
---8<--- "default/tutorials/tutorial1/fullpipeline_assembly.yml"
```


## metaSPAdes

Run metaspades to showcase the possibility to run/include other tools?
 
