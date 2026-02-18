## Creating grids

PineAPFEL can create PineAPPL interpolation grids pre-filled with analytically computed
coefficient functions from APFEL++. This is useful for structure function computations
where the hard-scattering kernels are known analytically and only the non-perturbative
input (PDFs or Fragmentation Functions) needs to be provided at convolution time.

The grid creation workflow requires three YAML configuration files:

1. A **grid card** that defines the process, observable, binning, channels, and perturbative orders
2. A **theory card** that specifies the QCD parameters (coupling, thresholds, perturbative order)
3. An **operator card** that defines the x-space interpolation grid and tabulation parameters

### Supported processes and observables

The `build_grid()` function currently supports the following combinations:

| Process | Observable | Current | APFEL++ initializer | Mass scheme |
|---------|-----------|---------|---------------------|-------------|
| DIS | F2 | NC | `InitializeF2NCObjectsZM` | Zero-mass |
| DIS | FL | NC | `InitializeFLNCObjectsZM` | Zero-mass |
| DIS | F3 | NC | `InitializeF3NCObjectsZM` | Zero-mass |
| SIA | F2 | NC | `InitializeF2NCObjectsZMT` | Zero-mass |
| SIA | FL | NC | `InitializeFLNCObjectsZMT` | Zero-mass |
| SIA | F3 | NC | `InitializeF3NCObjectsZMT` | Zero-mass |

All six combinations use the **zero-mass variable-flavour-number scheme** (ZM-VFNS), where
quarks are treated as massless above their respective thresholds.

!!! info "What is not yet supported"
    - **Charged-current (CC)** processes are not yet implemented. Setting `Current: CC`
      in the grid card will raise an error.
    - **SIDIS** (semi-inclusive DIS) grid filling is not supported. SIDIS grids can still
      be created with `create_grid()` and filled manually, but `build_grid()` does not
      handle the two-convolution structure.
    - **Massive coefficient functions** (FONLL, ACOT, S-ACOT) are not available.
    - **Polarised structure functions** (\(g_1\), \(g_4\), \(g_L\)) are not exposed
      through `build_grid()`, although APFEL++ does provide initializers for them.
    - **QCD+QED** coefficient functions are not implemented in grid filling (QED corrections
      are only available in the evolution step).

### Convolution types

Each process type requires a specific convolution type, which determines the type of
non-perturbative input the grid will be convoluted with:

| Process | `ConvolutionTypes` | Description |
|---------|-------------------|-------------|
| DIS | `[UNPOL_PDF]` | Unpolarised parton distribution functions |
| SIA | `[UNPOL_FF]` | Unpolarised fragmentation functions |
| SIDIS | `[UNPOL_PDF, UNPOL_FF]` | PDF for the initial state, FF for the final state |

Additionally, PineAPFEL's evolution module supports `POL_PDF` (longitudinally polarised PDFs)
for evolving grids that were filled externally with polarised coefficient functions.

### Perturbative orders

The coefficient functions are available at the following perturbative orders:

| `alpha_s` power | Label | F2/FL content | F3 content |
|:-:|:--:|---|---|
| 0 | LO | \(\delta(1-x)\) for quarks, 0 for gluon | \(\delta(1-x)\) for quarks, 0 for gluon |
| 1 | NLO | \(C_{2,\mathrm{NS}}^{(1)}\), \(C_{2,g}^{(1)}\) | \(C_{3,\mathrm{NS}}^{(1)}\), no gluon |
| 2 | NNLO | \(C_{2,\mathrm{NS}}^{(2)}\), \(C_{2,\mathrm{PS}}^{(2)}\), \(C_{2,g}^{(2)}\) | \(C_{3,\mathrm{NS}}^{(2)}\), no gluon |

The orders are specified in the grid card via the `Orders` field. Each order entry is a
5-element array `[alpha_s, alpha, log_xir, log_xif, log_xia]`. For pure QCD coefficient
functions, only the `alpha_s` power matters; the remaining entries should be set to 0.

Each entry stores the coefficient function at that **specific** power of \(\alpha_s\),
not the cumulative sum. For a complete NNLO prediction, all three orders (LO, NLO, NNLO)
must be listed so that the grid contains separate subgrids for each perturbative contribution.

!!! warning
    Orders beyond NNLO (`alpha_s > 2`) are silently skipped during grid filling, even
    though APFEL++ provides N3LO coefficient functions for some observables. This is a
    current limitation that may be lifted in a future release.

### Channel decomposition

The coefficient functions are decomposed in the **physical (PDG) basis**. Each channel in
the grid card maps to a combination of parton-level coefficient functions according to the
structure function decomposition.

For **F2** and **FL** (neutral current), the structure function is:

\[
F(x, Q^2) = \sum_{q} \mathcal{C}_q \otimes (q + \bar{q}) \;+\; \mathcal{C}_g \otimes g
\]

The per-channel coefficient functions are constructed from the APFEL++ operators
\(C_\mathrm{NS}\), \(C_\mathrm{S}\), and \(C_\mathrm{G}\) as follows:

| Channel | PIDs | Coefficient function |
|---------|------|---------------------|
| Quark \(q + \bar{q}\) | `[[q], [-q]]` | \(\mathcal{C}_q = e_q^2 \, C_\mathrm{NS} + \frac{\Sigma_\mathrm{ch}}{6}\,(C_\mathrm{S} - C_\mathrm{NS})\) |
| Gluon | `[[21]]` | \(\mathcal{C}_g = \Sigma_\mathrm{ch} \, C_\mathrm{G}\) |

where \(\Sigma_\mathrm{ch} = \sum_{i=1}^{n_f} e_i^2\) is the sum of electroweak charges
for the \(n_f\) active quark flavours at the given \(Q^2\), and the factor of 6 matches
the internal normalisation convention used in APFEL++'s `DISNCBasis`.

For **F3** (neutral current, parity-violating), only the non-singlet contributes:

| Channel | PIDs | Coefficient function |
|---------|------|---------------------|
| Quark \(q + \bar{q}\) | `[[q], [-q]]` | \(\mathcal{C}_q = e_q^2 \, C_\mathrm{NS}\) |
| Gluon | `[[21]]` | \(\mathcal{C}_g = 0\) |

!!! note
    F3 is a valence-type structure function. When using F3, the grid channels should use
    factors `[1.0, -1.0]` for quark channels (`[[q], [-q]]`) to produce the combination
    \(q - \bar{q}\), rather than `[1.0, 1.0]` which gives \(q + \bar{q}\).

---

### The grid card

The grid card defines the structure of the PineAPPL grid. Below is a complete reference
with all fields.

```yaml
# The scattering process.
# Supported values: DIS, SIA, SIDIS
Process: DIS

# The structure function observable.
# Supported values: F2, FL, F3
# Optional, defaults to F2.
Observable: F2

# The electroweak current type.
# Supported values: NC
# Optional, defaults to NC. CC is not yet supported.
Current: NC

# PID basis for the channel definitions.
# Supported values: PDG, EVOL
PidBasis: PDG

# Hadron PIDs involved in the process (one per convolution).
#   DIS:   [2212]        (proton)
#   SIA:   [211]         (pion)
#   SIDIS: [2212, 211]   (proton + pion)
HadronPids: [2212]

# Convolution types (one per convolution).
#   UNPOL_PDF  — unpolarised parton distributions
#   POL_PDF    — longitudinally polarised parton distributions
#   UNPOL_FF   — unpolarised fragmentation functions
#   POL_FF     — polarised fragmentation functions
ConvolutionTypes: [UNPOL_PDF]

# Perturbative orders as 5-element arrays:
# [alpha_s, alpha, log(xi_R), log(xi_F), log(xi_A)]
# For QCD-only coefficient functions, set all but alpha_s to 0.
Orders:
  - [0, 0, 0, 0, 0]   # O(alpha_s^0) = LO
  - [1, 0, 0, 0, 0]   # O(alpha_s^1) = NLO
  - [2, 0, 0, 0, 0]   # O(alpha_s^2) = NNLO

# Partonic channels. Each channel is a list of PID combinations
# with associated numerical factors.
# For DIS/SIA with one convolution, each combination has one PID.
Channels:
  - pids: [[2], [-2]]     # u + ubar
    factors: [1.0, 1.0]
  - pids: [[1], [-1]]     # d + dbar
    factors: [1.0, 1.0]
  - pids: [[21]]          # gluon
    factors: [1.0]

# Kinematic bins. Each bin is defined by lower and upper edges
# in each dimension:
#   DIS:   [Q^2, x]        (2 dimensions)
#   SIA:   [Q^2, z]        (2 dimensions)
#   SIDIS: [Q^2, x, z]     (3 dimensions)
Bins:
  - lower: [10.0, 0.001]
    upper: [100.0, 0.01]
  - lower: [100.0, 0.01]
    upper: [1000.0, 0.1]

# Bin normalisation factors (one per bin).
Normalizations: [1.0, 1.0]
```

#### DIS example

A DIS \(F_2\) grid up to NNLO with two \((Q^2, x)\) bins and three channels
(u+ubar, d+dbar, gluon):

```yaml
Process: DIS
Observable: F2
Current: NC
PidBasis: PDG
HadronPids: [2212]
ConvolutionTypes: [UNPOL_PDF]

Orders:
  - [0, 0, 0, 0, 0]
  - [1, 0, 0, 0, 0]
  - [2, 0, 0, 0, 0]

Channels:
  - pids: [[2], [-2]]
    factors: [1.0, 1.0]
  - pids: [[1], [-1]]
    factors: [1.0, 1.0]
  - pids: [[21]]
    factors: [1.0]

Bins:
  - lower: [10.0, 0.001]
    upper: [100.0, 0.01]
  - lower: [100.0, 0.01]
    upper: [1000.0, 0.1]

Normalizations: [1.0, 1.0]
```

#### SIA example

An SIA \(F_2\) grid for pion production up to NNLO. Note that the second kinematic
dimension is the hadron momentum fraction \(z\) instead of Bjorken \(x\),
and the convolution type is `UNPOL_FF`:

```yaml
Process: SIA
Observable: F2
Current: NC
PidBasis: PDG
HadronPids: [211]
ConvolutionTypes: [UNPOL_FF]

Orders:
  - [0, 0, 0, 0, 0]
  - [1, 0, 0, 0, 0]
  - [2, 0, 0, 0, 0]

Channels:
  - pids: [[2], [-2]]
    factors: [1.0, 1.0]
  - pids: [[1], [-1]]
    factors: [1.0, 1.0]
  - pids: [[21]]
    factors: [1.0]

Bins:
  - lower: [10.0, 0.2]
    upper: [100.0, 0.4]
  - lower: [100.0, 0.4]
    upper: [1000.0, 0.6]

Normalizations: [1.0, 1.0]
```

---

### Using the CLI

The `pineapfel-build` executable creates a filled PineAPPL grid from three YAML cards:

```text
pineapfel-build <grid.yaml> <theory.yaml> <operator.yaml> [-o output.pineappl.lz4]
```

If `-o` is not specified, the output filename is derived from the grid card by replacing
`.yaml` with `.pineappl.lz4`.

```bash
# Build a DIS F2 grid
pineapfel-build runcards/grid_dis.yaml runcards/theory.yaml runcards/operator.yaml

# Build an SIA grid with a custom output name
pineapfel-build runcards/grid_sia.yaml runcards/theory.yaml runcards/operator.yaml \
    -o sia_f2.pineappl.lz4
```

### Using the library

The same functionality is available programmatically through the `build_grid()` function:

```cpp
#include <pineapfel.h>
#include <pineappl_capi.h>
#include <iostream>

int main() {
    // 1. Load all three cards
    auto grid_def = pineapfel::load_grid_def("runcards/grid_dis.yaml");
    auto theory   = pineapfel::load_theory_card("runcards/theory.yaml");
    auto op_card  = pineapfel::load_operator_card("runcards/operator.yaml");

    // 2. Build and fill the grid
    auto* grid = pineapfel::build_grid(grid_def, theory, op_card);

    // 3. Write the grid
    pineappl_grid_write(grid, "dis_f2.pineappl.lz4");

    // 4. Cleanup
    pineappl_grid_delete(grid);
    return 0;
}
```

The returned grid can also be passed directly to `pineapfel::evolve()` to produce an
FK table in the same program:

```cpp
auto* grid    = pineapfel::build_grid(grid_def, theory, op_card);
auto* fktable = pineapfel::evolve(grid, theory, op_card);
pineappl_grid_write(fktable, "dis_f2.fk.pineappl.lz4");
pineappl_grid_delete(fktable);
pineappl_grid_delete(grid);
```

### Grid structure internals

Understanding the internal layout of the generated grid is useful for debugging and
for writing custom grid-filling code.

#### Node selection

The grid nodes are chosen automatically:

- **x/z nodes**: Taken from the APFEL++ joint interpolation grid, which is built from
  the `xgrid` definition in the operator card. This ensures consistency between the
  coefficient function grid and any subsequent evolution step.

- **Q^2 nodes**: Derived from the bin edges in the grid card. The bin boundaries in the
  first dimension (Q^2) are collected, and geometrically spaced intermediate points are
  added within each bin. Additionally, quark threshold values (\(m_q^2\)) that fall within
  the overall Q^2 range are included to capture flavour-threshold effects.

#### Subgrid layout

Each subgrid (one per combination of bin, perturbative order, and channel) is a
two-dimensional array of shape `[n_Q2, n_x]`, stored in row-major order. The
`node_values` vector concatenates the Q^2 nodes followed by the x/z nodes:

```
node_values = [Q^2_0, Q^2_1, ..., Q^2_{nq-1}, x_0, x_1, ..., x_{nx-1}]
```

For each Q^2 node, the coefficient function operator is evaluated at the bin's
x/z centre (geometric mean of the bin edges) to produce a distribution on the
APFEL++ joint grid, which populates one row of the subgrid.

### Programmatic grid definition

In addition to loading from YAML, you can construct a `GridDef` programmatically:

```cpp
pineapfel::GridDef def;
def.process    = pineapfel::ProcessType::DIS;
def.observable = pineapfel::Observable::F2;
def.current    = pineapfel::Current::NC;
def.pid_basis  = PINEAPPL_PID_BASIS_PDG;
def.hadron_pids        = {2212};
def.convolution_types  = {PINEAPPL_CONV_TYPE_UNPOL_PDF};

def.orders   = {{0, 0, 0, 0, 0}, {1, 0, 0, 0, 0}, {2, 0, 0, 0, 0}};  // LO + NLO + NNLO
def.channels = {
    {{{2}, {-2}},  {1.0, 1.0}},   // u + ubar
    {{{1}, {-1}},  {1.0, 1.0}},   // d + dbar
    {{{21}},       {1.0}},        // gluon
};
def.bins = {
    {{10.0, 0.001}, {100.0, 0.01}},
    {{100.0, 0.01}, {1000.0, 0.1}},
};
def.normalizations = {1.0, 1.0};

auto* grid = pineapfel::build_grid(def, theory, op_card);
```
