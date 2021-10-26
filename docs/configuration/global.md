# Global parameter settings

## Optional configuration of computational resources used for pipeline runs

The toolkit uses the following machine types (flavours) for running tools. All flavours can be optionally
adjusted by modifying the cpus and memory parameters. If for example the largest flavour is not available
in the infrastructure, `cpus` and `memory` parameters can be modified to fit the medium flavour. If larger
flavours are available, it makes especially sense to increase the `cpus` and `memory` values of the `large`
flavour to speed up for example assembly and read mapping.

Example Configuration:

```
resources:
  large:
    cpus: 28
    memory: 265
  medium:
    cpus: 14
    memory: 128
  small:
    cpus: 7
    memory: 16
  tiny:
    cpus: 1
    memory: 2
```

Additional flavours can be defined that can be used by methods that dynamically compute resources on tool error (see assembly module section).

Example:

```
resources:
  xlarge:
    cpus: 56
    memory: 512
  large:
    cpus: 28
    memory: 256
  medium:
    cpus: 14
    memory: 128
  small:
    cpus: 7
    memory: 16
  tiny:
    cpus: 1
    memory: 2
```
