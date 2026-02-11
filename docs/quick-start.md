`PineAPFEL` provides two main ways to evolve `PineAPPL` grids: it can be used either as a C++
library (see [the following section](library.md)) or through its command line interface (CLI).
In this section we describe how to perform the evolution using the CLI.

The executable `pineapfel-evolve` takes as input a `PineAPPL` grid file (typically with extension
`.pineappl.lz4`), together with two configuration files: a theory configuration file (`theory.yaml`)
and an operator configuration file (`operator.yaml`). The theory file specifies the physics setup
used for the evolution, such as the perturbative order, the treatment of the strong coupling, heavy-quark
thresholds, and scale choices, while the operator file defines the details of the evolution operator,
including the evolution basis, interpolation settings.

The program can be run as:
```bash
./build/pineapfel-evolve /path/to/grid.pineappl.lz4 runcards/theory.yaml runcards/operator.yaml
```

Upon successful execution, an FK table is produced in the same directory as the input grid. This FK
table contains the evolved grid and can subsequently be convoluted with parton distribution functions
for fast cross-section evaluations.
