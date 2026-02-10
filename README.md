# pineapfel

A C++ library and CLI tool for evolving [PineAPPL](https://github.com/NNPDF/pineappl)
interpolation grids into FK tables using [APFEL++](https://github.com/vbertone/apfelxx)
DGLAP evolution.

All physics parameters (coupling constants, thresholds, perturbative order, etc.)
are specified through two YAML configuration files — a **theory card** and an
**operator card** — rather than being hardcoded.

## Table of Contents

- [Dependencies](#dependencies)
- [Building](#building)
- [Quick Start](#quick-start)
- [Writing Configuration Cards](#writing-configuration-cards)
  - [Theory Card](#theory-card)
  - [Operator Card](#operator-card)
- [Using the CLI](#using-the-cli)
- [Using the Library](#using-the-library)
- [QCD vs QCD+QED Evolution](#qcd-vs-qcdqed-evolution)
- [File Structure](#file-structure)

## Dependencies

| Library | Purpose | Discovery |
|---------|---------|-----------|
| [APFEL++](https://github.com/vbertone/apfelxx) | DGLAP evolution kernels | `apfelxx-config` |
| [PineAPPL C API](https://github.com/NNPDF/pineappl) | Grid I/O and evolution interface | `pkg-config` |
| [yaml-cpp](https://github.com/jbeder/yaml-cpp) | YAML card parsing | `pkg-config` |
| CMake >= 3.14 | Build system | — |
| C++20 compiler | — | — |

## Building

```bash
cd pineapfel
mkdir build && cd build
cmake ..
make -j$(nproc)
```

This produces:
- `libpineapfel.a` — the static library
- `pineapfel-evolve` — the CLI executable

To install system-wide:

```bash
make install
```

> **Note:** If APFEL++ is installed in a non-standard location (e.g. `/usr/local/lib`),
> you may need to set `LD_LIBRARY_PATH` at runtime:
> ```bash
> export LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH
> ```

## Quick Start

Evolve a PineAPPL grid using the example configuration files:

```bash
./build/pineapfel-evolve /path/to/grid.pineappl.lz4 examples/theory.yaml examples/operator.yaml
```

This writes the FK table to `grid.fk.pineappl.lz4` (same directory as the input grid).

## Writing Configuration Cards

All physics and numerical parameters are split between two YAML files.

### Theory Card

The theory card specifies the physics of the DGLAP evolution. Here is a complete example with all
fields explained:

```yaml
# Initial evolution scale in GeV.
mu0: 1.0

# QCD perturbative order:
#   0 = LO
#   1 = NLO
#   2 = NNLO
PerturbativeOrder: 2

# Reference scale for alpha_s in GeV.
QRef: 91.1876

# Value of alpha_s(Q_ref)
AlphaQCDRef: 0.118

# Quark mass thresholds in GeV, one per active flavor.
# Ordering: up, down, strange, charm, bottom.
# Use 0.0 for massless quarks. The length determines the maximum number
# of active flavors (nf_max).
QuarkThresholds: [0.0, 0.0, 0.0, 1.4, 4.75]

# Output flavor PIDs written to the FK table.
# These define the pids_out dimension of the evolution operator.
# Standard PDG codes: 1=d, 2=u, 3=s, 4=c, 5=b, -1=dbar, ..., 21=gluon, 22=photon
Flavors: [-5, -4, -3, -2, -1, 1, 2, 3, 4, 5, 21]

# Set to true to use combined QCD+QED evolution.
# When false (default), pure QCD evolution is used.
QED: false

# --- QED-specific fields (ignored when QED: false) ---

# Value of alpha_em at the reference scale
AlphaQEDRef: 0.00729927

# Lepton mass thresholds in GeV.
# Leave empty to forbid lepton creation during evolution.
LeptonThresholds: []
```

#### Common configurations

**NNLO QCD with 5 flavors:**

```yaml
mu0: 1.0
PerturbativeOrder: 2
QRef: 91.1876
AlphaQCDRef: 0.118
QuarkThresholds: [0.0, 0.0, 0.0, 1.51, 4.92]
Flavors: [-5, -4, -3, -2, -1, 1, 2, 3, 4, 5, 21]
QED: false
```

**NLO QCD+QED with photon in the output basis:**

```yaml
mu0: 1.0
PerturbativeOrder: 1
QRef: 91.1876
AlphaQCDRef: 0.118
QuarkThresholds: [0.0, 0.0, 0.0, 1.51, 4.92]
Flavors: [-5, -4, -3, -2, -1, 1, 2, 3, 4, 5, 21, 22]
QED: true
AlphaQEDRef: 0.00729927
LeptonThresholds: []
```

### Operator Card

The operator card controls numerical aspects of the evolution: the x-space interpolation grid,
the tabulation of running couplings, and scale-variation factors.

```yaml
# x-space interpolation grid.
# Each sub-grid is defined by:
#   n_knots:     number of interpolation knots
#   x_min:       lower edge of the sub-grid in Bjorken x
#   poly_degree: polynomial interpolation degree
#
# The sub-grids are joined into a single APFEL++ Grid object.
# Add more sub-grids for finer resolution in specific x regions.
xgrid:
  - n_knots: 80
    x_min: 1.0e-5
    poly_degree: 3
  - n_knots: 40
    x_min: 1.0e-1
    poly_degree: 3

# Tabulation parameters for the running coupling(s).
# These are passed to the apfel::TabulateObject constructor:
#   n_points:     number of interpolation points
#   q_min:        minimum scale in GeV for the tabulation
#   n_steps:      number of steps in the tabulation grid
#   interp_degree: interpolation degree for the tabulated object
tabulation:
  n_points: 200
  q_min: 0.9
  n_steps: 13001
  interp_degree: 3

# Scale variation factors: [xi_R, xi_F, xi_A]
#   xi_R: renormalization scale ratio (mu_R / Q)
#   xi_F: factorization scale ratio  (mu_F / Q)
#   xi_A: fragmentation scale ratio  (mu_A / Q)
# Use [1.0, 1.0, 1.0] for central predictions.
xi: [1.0, 1.0, 1.0]
```

## Using the CLI

```
pineapfel-evolve <grid.pineappl.lz4> <theory.yaml> <operator.yaml> [-o output.pineappl.lz4]
```

| Argument | Description |
|----------|-------------|
| `grid.pineappl.lz4` | Input PineAPPL interpolation grid |
| `theory.yaml` | Theory card (physics parameters) |
| `operator.yaml` | Operator card (numerical parameters) |
| `-o output` | *(Optional)* Output path. Default: `<grid-name>.fk.pineappl.lz4` |

**Examples:**

```bash
# Default output name
pineapfel-evolve grid.pineappl.lz4 theory.yaml operator.yaml

# Custom output name
pineapfel-evolve grid.pineappl.lz4 theory.yaml operator.yaml -o my_fktable.pineappl.lz4
```

## Using the Library

Include the single convenience header to access the full API:

```cpp
#include <pineapfel/pineapfel.h>
```

### Minimal example

```cpp
#include <pineapfel/pineapfel.h>
#include <pineappl_capi.h>
#include <iostream>

int main() {
    // 1. Load configuration cards from YAML
    auto theory  = pineapfel::load_theory_card("theory.yaml");
    auto op_card = pineapfel::load_operator_card("operator.yaml");

    // 2. Read the PineAPPL grid
    auto* grid = pineappl_grid_read("grid.pineappl.lz4");

    // 3. Evolve into an FK table
    auto* fktable = pineapfel::evolve(grid, theory, op_card);

    // 4. Write out the FK table
    pineappl_grid_write(fktable, "fktable.pineappl.lz4");

    // 5. Cleanup
    pineappl_grid_delete(fktable);
    pineappl_grid_delete(grid);

    return 0;
}
```

### Building your own program against `pineapfel`

If `pineapfel` is installed, you can link against it in your `CMakeLists.txt`:

```cmake
find_library(PINEAPFEL_LIB pineapfel)
target_link_libraries(my_program PRIVATE ${PINEAPFEL_LIB})
target_include_directories(my_program PRIVATE /path/to/pineapfel/include)
```

Or compile directly with the source tree:

```cmake
add_subdirectory(pineapfel)
target_link_libraries(my_program PRIVATE pineapfel)
```

### Constructing cards programmatically

You do not need YAML files — the structs can be filled directly in C++:

```cpp
#include <pineapfel/pineapfel.h>
#include <pineappl_capi.h>

int main() {
    // Build a theory card in code
    pineapfel::TheoryCard theory;
    theory.mu0            = 1.0;
    theory.pert_ord       = 2;            // NNLO
    theory.q_ref          = 91.1876;      // M_Z
    theory.alpha_qcd_ref  = 0.118;
    theory.quark_thresholds = {0.0, 0.0, 0.0, 1.51, 4.92};
    theory.flavors        = {-5, -4, -3, -2, -1, 1, 2, 3, 4, 5, 21};
    theory.qed            = false;

    // Build an operator card in code
    pineapfel::OperatorCard op_card;
    op_card.xgrid = {
        {80, 1e-5, 3},   // {n_knots, x_min, poly_degree}
        {40, 1e-1, 3},
    };
    op_card.tabulation = {200, 0.9, 13001, 3};  // {n_points, q_min, n_steps, interp_degree}
    op_card.xi = {1.0, 1.0, 1.0};

    // Load grid and evolve
    auto* grid = pineappl_grid_read("grid.pineappl.lz4");
    auto* fktable = pineapfel::evolve(grid, theory, op_card);

    pineappl_grid_write(fktable, "fktable.pineappl.lz4");
    pineappl_grid_delete(fktable);
    pineappl_grid_delete(grid);

    return 0;
}
```

## QCD vs QCD+QED Evolution

The evolution mode is controlled by the `QED` field in the theory card:

| Mode | `QED` | Basis size | DGLAP objects | Coupling |
|------|-------|-----------|---------------|----------|
| Pure QCD | `false` | 13 x 13 | `apfel::DglapObjects` | `apfel::AlphaQCD` |
| QCD+QED | `true` | 20 x 20 | `apfel::DglapObjectsQCDQED` | `apfel::AlphaQCDQED` |

**QCD mode** supports three convolution types (determined automatically from the grid):
- `UNPOL_PDF` — unpolarised PDFs
- `POL_PDF` — longitudinally polarised PDFs
- `UNPOL_FF` — unpolarised fragmentation functions

**QCD+QED mode** currently supports:
- `UNPOL_PDF` — unpolarised PDFs with QED corrections

The library reads the convolution types from the grid metadata and initializes the appropriate
DGLAP objects automatically.
