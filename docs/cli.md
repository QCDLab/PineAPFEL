## Using the CLI

PineAPFEL provides two command-line executables:

- `pineapfel-build` — creates a PineAPPL grid filled with APFEL++ coefficient functions
- `pineapfel-evolve` — evolves a PineAPPL grid into an FK table

### `pineapfel-build`

Creates a PineAPPL grid from a grid card, theory card, and operator card. The grid is
filled with analytically computed coefficient functions for the specified process and
observable (see [Grid creation](grid-creation.md) for details on supported processes).

```text
pineapfel-build <grid.yaml> <theory.yaml> <operator.yaml> [-o output.pineappl.lz4]
```

If `-o` is not specified, the output filename is derived from the grid card by replacing
`.yaml` with `.pineappl.lz4`.

#### Examples

```bash
# Build a DIS F2 grid (output: runcards/grid_dis.pineappl.lz4)
pineapfel-build runcards/grid_dis.yaml runcards/theory.yaml runcards/operator.yaml

# Build an SIA grid with a custom output name
pineapfel-build runcards/grid_sia.yaml runcards/theory.yaml runcards/operator.yaml \
    -o sia_f2.pineappl.lz4

# Build a SIDIS F2 grid (two-convolution: PDF ⊗ FF)
pineapfel-build runcards/grid_sidis.yaml runcards/theory.yaml runcards/operator.yaml \
    -o sidis_f2.pineappl.lz4

# Build a DIS F2 grid in the FFN massive scheme
# (requires MassScheme: FFN in the grid card and HeavyQuarkMasses in the theory card)
pineapfel-build runcards/grid_dis_ffn.yaml runcards/theory.yaml runcards/operator.yaml \
    -o dis_f2_ffn.pineappl.lz4

# Build a DIS F2 FONLL grid (F_ZM + F_FFN)
pineapfel-build runcards/grid_dis_fonll.yaml runcards/theory.yaml runcards/operator.yaml \
    -o dis_f2_fonll.pineappl.lz4
```

### `pineapfel-evolve`

Evolves an existing PineAPPL interpolation grid into an FK table using APFEL++ DGLAP
evolution operators. The input grid can be one produced by `pineapfel-build` or any
other valid PineAPPL grid.

```text
pineapfel-evolve <grid.pineappl.lz4> <theory.yaml> <operator.yaml> [-o output.pineappl.lz4]
```

If `-o` is not specified, the output is written to `<input>.fk.pineappl.lz4`.

#### Examples

```bash
# Default output name
pineapfel-evolve grid.pineappl.lz4 theory.yaml operator.yaml

# Custom output name
pineapfel-evolve grid.pineappl.lz4 theory.yaml operator.yaml -o my_fktable.pineappl.lz4
```

### Full pipeline

A typical workflow chains grid creation and evolution. The same commands work for DIS,
SIA, and SIDIS — only the grid card differs:

```bash
# DIS F2: create and evolve
pineapfel-build runcards/grid_dis.yaml runcards/theory.yaml runcards/operator.yaml \
    -o dis_f2.pineappl.lz4
pineapfel-evolve dis_f2.pineappl.lz4 runcards/theory.yaml runcards/operator.yaml \
    -o dis_f2.fk.pineappl.lz4

# SIDIS F2: create and evolve (two-convolution grid)
pineapfel-build runcards/grid_sidis.yaml runcards/theory.yaml runcards/operator.yaml \
    -o sidis_f2.pineappl.lz4
pineapfel-evolve sidis_f2.pineappl.lz4 runcards/theory.yaml runcards/operator.yaml \
    -o sidis_f2.fk.pineappl.lz4
```

The resulting FK table can be convoluted with PDFs (and FFs for SIDIS) using the
`pineappl` CLI or any PineAPPL-compatible tool.
