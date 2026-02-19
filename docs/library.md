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
