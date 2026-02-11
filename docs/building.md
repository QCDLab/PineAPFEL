The following provides detailed instructions for building and installing `PineAPFEL`
and its APIs and tools. `PineAPFEL` depends on the following libraries:

- [PineAPPL](https://github.com/NNPDF/pineappl) for producing and reading fast interpolation
    grids
- [APFEL++](https://github.com/vbertone/apfelxx) for the DGLAP evolution kernels
- [YAML-CPP](https://github.com/jbeder/yaml-cpp) for parsing the theory and operator YAML
    cards

### CMake

```bash
cd pineapfel
mkdir build && cd build
cmake ..
make -j"$(nproc)"
```

To install system-wide:

```bash
make install
```

### Meson

```bash
cd pineapfel
meson setup build-meson
meson compile -C build-meson
```

To install system-wide:

```bash
meson install -C build-meson
```

!!! note
    If APFEL++ is installed in a non-standard location, you may need to set `LD_LIBRARY_PATH`
    at runtime:

    ```bash
    export LD_LIBRARY_PATH=/path/to/lib:$LD_LIBRARY_PATH
    ```

