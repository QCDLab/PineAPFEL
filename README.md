<h1 align="center">PineAPFEL</h1>

<p align="justify">
  A C++ library and CLI tool for evolving <a href="https://github.com/NNPDF/pineappl">PineAPPL</a>
  interpolation grids into FK tables using <a href="https://github.com/vbertone/apfelxx">APFEL++</a>
  DGLAP evolution.
</p>

<p align="justify">
  All physics parameters (coupling constants, thresholds, perturbative order, etc.)
  are specified through two YAML configuration files — a <b>theory card</b> and an
  <b>operator card</b> — rather than being hardcoded.
</p>

<h2>Table of Contents</h2>

<ul>
  <li><a href="#dependencies">Dependencies</a></li>
  <li><a href="#building">Building</a></li>
  <li><a href="#quick-start">Quick Start</a></li>
  <li><a href="#writing-configuration-cards">Writing Configuration Cards</a>
    <ul>
      <li><a href="#theory-card">Theory Card</a></li>
      <li><a href="#operator-card">Operator Card</a></li>
    </ul>
  </li>
  <li><a href="#using-the-cli">Using the CLI</a></li>
  <li><a href="#using-the-library">Using the Library</a></li>
  <li><a href="#qcd-vs-qcdqed-evolution">QCD vs QCD+QED Evolution</a></li>
  <li><a href="#file-structure">File Structure</a></li>
</ul>

<h2>Dependencies</h2>

<div align="center">
<table>
  <tr>
    <th>Library</th>
    <th>Purpose</th>
    <th>Discovery</th>
  </tr>
  <tr>
    <td><a href="https://github.com/vbertone/apfelxx">APFEL++</a></td>
    <td>DGLAP evolution kernels</td>
    <td><code>apfelxx-config</code></td>
  </tr>
  <tr>
    <td><a href="https://github.com/NNPDF/pineappl">PineAPPL C API</a></td>
    <td>Grid I/O and evolution interface</td>
    <td><code>pkg-config</code></td>
  </tr>
  <tr>
    <td><a href="https://github.com/jbeder/yaml-cpp">yaml-cpp</a></td>
    <td>YAML card parsing</td>
    <td><code>pkg-config</code></td>
  </tr>
  <tr>
    <td>CMake >= 3.14</td>
    <td>Build system</td>
    <td>—</td>
  </tr>
  <tr>
    <td>C++20 compiler</td>
    <td>—</td>
    <td>—</td>
  </tr>
</table>
</div>

<h2>Building</h2>

```bash
cd pineapfel
mkdir build && cd build
cmake ..
make -j$(nproc)
```

<p align="justify">
  This produces a static library (<code>libpineapfel.a</code>) and a CLI executable (<code>pineapfel-evolve</code>).
</p>

<p align="justify">
  To install system-wide:
</p>

```bash
make install
```

<blockquote>
<p align="justify">
  <b>Note:</b> If APFEL++ is installed in a non-standard location, you may need to set <code>LD_LIBRARY_PATH</code> at
  runtime:
</p>

```bash
export LD_LIBRARY_PATH=/path/to/lib:$LD_LIBRARY_PATH
```
</blockquote>

<h2>Quick Start</h2>

<p align="justify">
  Evolve a PineAPPL grid using the example configuration files:
</p>

```bash
./build/pineapfel-evolve /path/to/grid.pineappl.lz4 examples/theory.yaml examples/operator.yaml
```

<p align="justify">
  This writes the FK table to <code>grid.fk.pineappl.lz4</code> (same directory as the input grid).
</p>

<h2>Writing Configuration Cards</h2>

<p align="justify">
  All physics and numerical parameters are split between two YAML files.
</p>

<h3>Theory Card</h3>

<p align="justify">
  The theory card specifies the physics of the DGLAP evolution. Here is a complete example with all
  fields explained:
</p>

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

<h4>Common configurations</h4>

<p align="justify">
  <b>NNLO QCD with 5 flavors:</b>
</p>

```yaml
mu0: 1.0
PerturbativeOrder: 2
QRef: 91.1876
AlphaQCDRef: 0.118
QuarkThresholds: [0.0, 0.0, 0.0, 1.51, 4.92]
Flavors: [-5, -4, -3, -2, -1, 1, 2, 3, 4, 5, 21]
QED: false
```

<p align="justify">
  <b>NLO QCD+QED with photon in the output basis:</b>
</p>

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

<h3>Operator Card</h3>

<p align="justify">
  The operator card controls numerical aspects of the evolution: the x-space interpolation grid,
  the tabulation of running couplings, and scale-variation factors.
</p>

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
#   xi_F: factorization scale ratio   (mu_F / Q)
#   xi_A: fragmentation scale ratio   (mu_A / Q)
# Use [1.0, 1.0, 1.0] for central predictions.
xi: [1.0, 1.0, 1.0]
```

<h2>Using the CLI</h2>

```
pineapfel-evolve <grid.pineappl.lz4> <theory.yaml> <operator.yaml> [-o output.pineappl.lz4]
```

<div align="center">
<table>
  <tr>
    <th>Argument</th>
    <th>Description</th>
  </tr>
  <tr>
    <td><code>grid.pineappl.lz4</code></td>
    <td>Input PineAPPL interpolation grid</td>
  </tr>
  <tr>
    <td><code>theory.yaml</code></td>
    <td>Theory card (physics parameters)</td>
  </tr>
  <tr>
    <td><code>operator.yaml</code></td>
    <td>Operator card (numerical parameters)</td>
  </tr>
  <tr>
    <td><code>-o output</code></td>
    <td>(Optional) Output path. Default: <code>&lt;grid-name&gt;.fk.pineappl.lz4</code></td>
  </tr>
</table>
</div>

<p align="justify">
  <b>Examples:</b>
</p>

```bash
# Default output name
pineapfel-evolve grid.pineappl.lz4 theory.yaml operator.yaml

# Custom output name
pineapfel-evolve grid.pineappl.lz4 theory.yaml operator.yaml -o my_fktable.pineappl.lz4
```

<h2>Using the Library</h2>

<p align="justify">
  Include the single convenience header to access the full API:
</p>

```cpp
#include <pineapfel/pineapfel.h>
```

<h3>Minimal example</h3>

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

<h3>Building your own program against <code>pineapfel</code></h3>

<p align="justify">
  If <code>pineapfel</code> is installed, you can link against it in your <code>CMakeLists.txt</code>:
</p>

```cmake
find_library(PINEAPFEL_LIB pineapfel)
target_link_libraries(my_program PRIVATE ${PINEAPFEL_LIB})
target_include_directories(my_program PRIVATE /path/to/pineapfel/include)
```

<p align="justify">
  Or compile directly with the source tree:
</p>

```cmake
add_subdirectory(pineapfel)
target_link_libraries(my_program PRIVATE pineapfel)
```

<h3>Constructing cards programmatically</h3>

<p align="justify">
  You do not need YAML files — the structs can be filled directly in C++:
</p>

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

<h2>QCD vs QCD+QED Evolution</h2>

<p align="justify">
  The evolution mode is controlled by the <code>QED</code> field in the theory card:
</p>

<div align="center">
<table>
  <tr>
    <th>Mode</th>
    <th><code>QED</code></th>
    <th>Basis size</th>
    <th>DGLAP objects</th>
    <th>Coupling</th>
  </tr>
  <tr>
    <td>Pure QCD</td>
    <td><code>false</code></td>
    <td>13 x 13</td>
    <td><code>apfel::DglapObjects</code></td>
    <td><code>apfel::AlphaQCD</code></td>
  </tr>
  <tr>
    <td>QCD+QED</td>
    <td><code>true</code></td>
    <td>20 x 20</td>
    <td><code>apfel::DglapObjectsQCDQED</code></td>
    <td><code>apfel::AlphaQCDQED</code></td>
  </tr>
</table>
</div>

<p align="justify">
  <b>QCD mode</b> supports three convolution types (determined automatically from the grid):
</p>

<ul>
  <li><code>UNPOL_PDF</code> — unpolarised PDFs</li>
  <li><code>POL_PDF</code> — longitudinally polarised PDFs</li>
  <li><code>UNPOL_FF</code> — unpolarised fragmentation functions</li>
</ul>

<p align="justify">
  <b>QCD+QED mode</b> currently supports:
</p>

<ul>
  <li><code>UNPOL_PDF</code> — unpolarised PDFs with QED corrections</li>
</ul>

<p align="justify">
  The library reads the convolution types from the grid metadata and initializes the appropriate
  DGLAP objects automatically.
</p>
