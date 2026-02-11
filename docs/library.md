## Using the library

Include the single convenience header to access the full API:

```cpp
#include <pineapfel.h>
```

### Minimal example

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

