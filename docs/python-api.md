## Python API

`PineAPFEL` ships a Python package — `pineapfel` — built with
[pybind11](https://github.com/pybind/pybind11) and compiled as part of the normal
Meson build. It exposes the same card-loading, grid-building, and evolution
functions as the C++ library, and returns grids that can be directly handed to
[PineAPPL](https://github.com/NNPDF/pineappl) for convolution and analysis.

See [Building](building.md#python-api) for installation instructions.

---

### Grid creation

Build a `PineAPPL` grid filled with `APFEL++` coefficient functions:

```python
import pineapfel

# 1. Load the three YAML cards
theory  = pineapfel.load_theory_card("runcards/theory.yaml")
op_card = pineapfel.load_operator_card("runcards/operator.yaml")
gdef    = pineapfel.load_grid_def("runcards/grid_dis.yaml")

# 2. Build and fill the grid with coefficient functions
grid = pineapfel.build_grid(gdef, theory, op_card)

# 3. Write to disk
grid.write("dis_f2.pineappl.lz4")
```

The same `build_grid()` call works for all supported processes; the relevant
flags are read from the grid card:

```python
# SIDIS (PDF ⊗ FF, two convolution functions)
sidis_def = pineapfel.load_grid_def("runcards/grid_sidis.yaml")
sidis_grid = pineapfel.build_grid(sidis_def, theory, op_card)
sidis_grid.write("sidis_f2.pineappl.lz4")

# Polarized DIS g₁
pol_def  = pineapfel.load_grid_def("runcards/grid_dis_pol.yaml")
pol_grid = pineapfel.build_grid(pol_def, theory, op_card)
pol_grid.write("dis_g1.pineappl.lz4")

# Fixed-Flavour Number (FFN) massive scheme
ffn_def  = pineapfel.load_grid_def("runcards/grid_dis_ffn.yaml")
ffn_grid = pineapfel.build_grid(ffn_def, theory, op_card)
ffn_grid.write("dis_f2_ffn.pineappl.lz4")

# FONLL combined scheme
fonll_def  = pineapfel.load_grid_def("runcards/grid_dis_fonll.yaml")
fonll_grid = pineapfel.build_grid(fonll_def, theory, op_card)
fonll_grid.write("dis_f2_fonll.pineappl.lz4")
```

---

### Fixing a convolution

`pineapfel.fix_convolution` folds one convolution slot of a grid into a PDF/FF
evaluation, returning a new grid with one fewer convolution. This mirrors the
`FixConvolutions` YAML key available in the CLI workflow, but accepts **any Python
callable** as the PDF evaluator — a NeoPDF wrapper, a `lhapdf`-python object, or a
custom function.

```python
import pineapfel

theory  = pineapfel.load_theory_card("runcards/theory.yaml")
op_card = pineapfel.load_operator_card("runcards/operator.yaml")
sidis_def = pineapfel.load_grid_def("runcards/grid_sidis.yaml")

# Build the two-convolution SIDIS grid (PDF ⊗ FF)
sidis_grid = pineapfel.build_grid(sidis_def, theory, op_card)

# Define an xfx callable — here using a NeoPDF set via its Python interface
from neopdf.pdf import PDF as NeoPDF # NeoPDF Python bindings
pdf = NeoPDF.mkPDF("NNPDF40_nnlo_as_01180", member=0)

def xfx(pid: int, x: float, q2: float) -> float:
    return pdf.xfxQ2(pid, x, q2)  # must return x * f(x, Q²)

# Fix convolution slot 0 (UNPOL_PDF); result has one convolution (UNPOL_FF)
ff_grid = pineapfel.fix_convolution(sidis_grid, conv_idx=0, xfx=xfx, xi=1.0)

# Write the FF-only grid; can be evolved like any single-convolution grid
ff_grid.write("sidis_ff_only.pineappl.lz4")
```

The function signature is:

```python
pineapfel.fix_convolution(
    grid,        # Grid — source grid (not modified)
    conv_idx,    # int  — 0-based convolution slot to fix
    xfx,         # callable(pid: int, x: float, q2: float) -> float
    xi = 1.0,    # float — factorization-scale variation factor (default 1.0)
) -> Grid
```

The `FixConvolutionDef` data class is also exposed for use with programmatic
`GridDef` construction:

```python
fc = pineapfel.FixConvolutionDef(
    index=0,
    pdf_set="NNPDF40_nnlo_as_01180",
    pdf_member=0,   # optional, default 0
    xi=1.0,         # optional, default 1.0
)

gdef = pineapfel.load_grid_def("runcards/grid_sidis.yaml")
gdef.fix_convolutions = [fc]
```

!!! note
    `build_grid()` does not apply `fix_convolutions` entries from the `GridDef` —
    that step is the responsibility of the caller (or of `pineapfel-build` in the CLI
    workflow). This keeps the build step pure and lets you inspect the full
    two-convolution grid before folding.

---

### Evolution

Evolve an existing grid into an FK table:

```python
import pineapfel

theory  = pineapfel.load_theory_card("runcards/theory.yaml")
op_card = pineapfel.load_operator_card("runcards/operator.yaml")

# Read an existing PineAPPL grid
grid = pineapfel.Grid.read("dis_f2.pineappl.lz4")

# Evolve into an FK table
fktable = pineapfel.evolve(grid, theory, op_card)
fktable.write("dis_f2.fk.pineappl.lz4")
```

---

### Convolution using PineAPPL API

The `Grid` object returned by `build_grid` and `evolve` can be written to disk
and read back by [PineAPPL](https://github.com/NNPDF/pineappl) for convolution
with PDFs:

```python
import lhapdf
import numpy as np
import pineapfel

from pineappl.grid import Grid
from pineappl.convolutions import Conv, ConvType

# ── Build the grid
theory  = pineapfel.load_theory_card("runcards/theory.yaml")
op_card = pineapfel.load_operator_card("runcards/operator.yaml")
gdef    = pineapfel.load_grid_def("runcards/grid_dis.yaml")
grid    = pineapfel.build_grid(gdef, theory, op_card)
grid.write("dis_f2.pineappl.lz4")

# ── Read back with PineAPPL
pg = Grid.read("dis_f2.pineappl.lz4")

# ── Set up the PDF set
pdf = lhapdf.mkPDF("NNPDF40_nnlo_as_01180", 0)

def xfx(pid, x, q2):
    return pdf.xfxQ2(pid, x, q2)

def alphas(q2):
    return pdf.alphasQ2(q2)

# ── Convolve
ct  = ConvType(polarized=False, time_like=False)
c   = Conv(ct, 2212)                             # proton (PDG id 2212)
res = pg.convolve(pdg_convs=[c], xfxs=[xfx], alphas=alphas)
print("F2 per bin:", res)
```

---

### Programmatic grid definition

Cards can also be constructed in Python without YAML files:

```python
import pineapfel

theory = pineapfel.TheoryCard()
theory.mu0           = 1.0
theory.pert_ord      = 2
theory.q_ref         = 91.2
theory.alpha_qcd_ref = 0.118
theory.quark_thresholds = [0.3, 1.5, 4.5, 100.0, 100.0]
theory.flavors          = [1, 2, 3, 4, 5, 6]
theory.ckm              = [1.0, 0.0, 0.0] * 3   # diagonal CKM
theory.qed              = False
theory.alpha_qed_ref    = 1.0 / 137.0
theory.lepton_thresholds = []
theory.heavy_quark_masses = [0.3, 1.5, 4.5, 100.0, 100.0, 100.0]

op_card = pineapfel.OperatorCard()
op_card.xgrid      = [pineapfel.SubGridDef(100, 1e-7, 3)]
op_card.tabulation = pineapfel.TabulationParams(50, 1.0, 200, 3)
op_card.xi         = [1.0, 1.0, 1.0]

gdef = pineapfel.GridDef()
gdef.process           = pineapfel.ProcessType.DIS
gdef.observable        = pineapfel.Observable.F2
gdef.current           = pineapfel.Current.NC
gdef.mass_scheme       = pineapfel.MassScheme.ZM
gdef.pid_basis         = pineapfel.PidBasis.PDG
gdef.hadron_pids       = [2212]
gdef.convolution_types = [pineapfel.ConvolutionType.UNPOL_PDF]
gdef.orders            = [pineapfel.OrderDef(0, 0, 0, 0, 0),
                          pineapfel.OrderDef(1, 0, 0, 0, 0)]
gdef.bins = [
    pineapfel.BinDef([10.0,  0.001], [100.0,  0.01]),
    pineapfel.BinDef([100.0, 0.01 ], [1000.0, 0.1 ]),
]
gdef.normalizations = [1.0, 1.0]

grid = pineapfel.build_grid(gdef, theory, op_card)
grid.write("my_dis_f2.pineappl.lz4")
```
