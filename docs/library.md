## Using the library

Include the single convenience header to access the full API:

```cpp
#include <pineapfel.h>
```

### API overview

PineAPFEL exposes three main modules:

| Module | Header | Functions |
|--------|--------|-----------|
| Cards | `cards.h` | `load_theory_card()`, `load_operator_card()` |
| Grid creation | `grid_gen.h`, `fill.h` | `load_grid_def()`, `derive_channels()`, `create_grid()`, `build_grid()` |
| Evolution | `evolution.h` | `evolve()` |
| SIDIS coefficients | `sidis_api.h` | `init_sidis()`, `SidisCoeffs` |

`sidis_api.h` is an internal header used by `build_grid()` to load the SIDIS coefficient
functions from APFEL++ without exposing `<apfel/SIDIS.h>` to other translation units.
You do not need to include it in application code.

### Grid creation example

Build a PineAPPL grid filled with APFEL++ coefficient functions:

```cpp
#include <pineapfel.h>
#include <pineappl_capi.h>
#include <iostream>

int main() {
    // 1. Load the three YAML cards
    auto grid_def = pineapfel::load_grid_def("grid_dis.yaml");
    auto theory   = pineapfel::load_theory_card("theory.yaml");
    auto op_card  = pineapfel::load_operator_card("operator.yaml");

    // 2. Build and fill the grid with coefficient functions
    auto* grid = pineapfel::build_grid(grid_def, theory, op_card);

    // 3. Write the grid
    pineappl_grid_write(grid, "dis_f2.pineappl.lz4");
    pineappl_grid_delete(grid);

    return 0;
}
```

The same `build_grid()` call works unchanged for SIDIS, polarized grids, and massive
schemes — the relevant flags are read from the grid card:

```cpp
// SIDIS (two-convolution: PDF ⊗ FF)
auto grid_def = pineapfel::load_grid_def("grid_sidis.yaml");  // Process: SIDIS
auto* grid    = pineapfel::build_grid(grid_def, theory, op_card);
pineappl_grid_write(grid, "sidis_f2.pineappl.lz4");
pineappl_grid_delete(grid);

// Polarized DIS g1 (Polarized: true in card)
auto pol_def = pineapfel::load_grid_def("grid_dis_pol.yaml");
auto* pol_grid = pineapfel::build_grid(pol_def, theory, op_card);
pineappl_grid_write(pol_grid, "dis_g1.pineappl.lz4");
pineappl_grid_delete(pol_grid);

// FFN DIS F2 (MassScheme: FFN in card, HeavyQuarkMasses in theory card)
auto ffn_def = pineapfel::load_grid_def("grid_dis_ffn.yaml");
auto* ffn_grid = pineapfel::build_grid(ffn_def, theory, op_card);
pineappl_grid_write(ffn_grid, "dis_f2_ffn.pineappl.lz4");
pineappl_grid_delete(ffn_grid);

// FONLL DIS F2 (MassScheme: FONLL in card)
auto fonll_def = pineapfel::load_grid_def("grid_dis_fonll.yaml");
auto* fonll_grid = pineapfel::build_grid(fonll_def, theory, op_card);
pineappl_grid_write(fonll_grid, "dis_f2_fonll.pineappl.lz4");
pineappl_grid_delete(fonll_grid);
```

When building programmatically, polarization is inferred from `convolution_types`:
setting the first entry to `PINEAPPL_CONV_TYPE_POL_PDF` selects polarized coefficient
functions. Use `grid_def.mass_scheme` to select a heavy-quark scheme:

```cpp
// Polarized DIS g1: POL_PDF in convolution_types drives the coefficient function choice
pineapfel::GridDef def;
def.process    = pineapfel::ProcessType::DIS;
def.observable = pineapfel::Observable::F2;  // interpreted as g1
def.current    = pineapfel::Current::NC;
def.pid_basis  = PINEAPPL_PID_BASIS_PDG;
def.hadron_pids       = {2212};
def.convolution_types = {PINEAPPL_CONV_TYPE_POL_PDF};
def.orders = {{0, 0, 0, 0, 0}, {1, 0, 0, 0, 0}};
def.bins = {{{10.0, 0.001}, {100.0, 0.01}}};
def.normalizations = {1.0};

auto* pol_grid = pineapfel::build_grid(def, theory, op_card);

// FFN DIS F2 (theory card must contain HeavyQuarkMasses with 6 entries)
pineapfel::GridDef ffn_def;
ffn_def.process      = pineapfel::ProcessType::DIS;
ffn_def.observable   = pineapfel::Observable::F2;
ffn_def.current      = pineapfel::Current::NC;
ffn_def.mass_scheme  = pineapfel::MassScheme::FFN;   // <-- key field
ffn_def.pid_basis    = PINEAPPL_PID_BASIS_PDG;
ffn_def.hadron_pids       = {2212};
ffn_def.convolution_types = {PINEAPPL_CONV_TYPE_UNPOL_PDF};
ffn_def.orders = {{0, 0, 0, 0, 0}, {1, 0, 0, 0, 0}, {2, 0, 0, 0, 0}};
ffn_def.bins = {{{10.0, 0.001}, {100.0, 0.01}}};
ffn_def.normalizations = {1.0};

auto* ffn_grid = pineapfel::build_grid(ffn_def, theory, op_card);
```

See [Grid creation](grid-creation.md) for details on the grid card format and
supported processes.

### Evolution example

Evolve an existing PineAPPL grid into an FK table:

```cpp
#include <pineapfel.h>
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

### Full pipeline example

Create a grid and evolve it in one program:

```cpp
#include <pineapfel.h>
#include <pineappl_capi.h>

int main() {
    auto grid_def = pineapfel::load_grid_def("grid_dis.yaml");
    auto theory   = pineapfel::load_theory_card("theory.yaml");
    auto op_card  = pineapfel::load_operator_card("operator.yaml");

    // Build the coefficient function grid
    auto* grid = pineapfel::build_grid(grid_def, theory, op_card);

    // Evolve into an FK table
    auto* fktable = pineapfel::evolve(grid, theory, op_card);

    // Write out
    pineappl_grid_write(fktable, "dis_f2.fk.pineappl.lz4");

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
