The following provides detailed instructions for building and installing `PineAPFEL`
and its APIs and tools. `PineAPFEL` depends on the following libraries:

- [PineAPPL](https://github.com/NNPDF/pineappl) for producing and reading fast interpolation
    grids
- [APFEL++](https://github.com/vbertone/apfelxx) for the DGLAP evolution kernels
- [YAML-CPP](https://github.com/jbeder/yaml-cpp) for parsing the theory and operator YAML
    cards

## C++ library and CLI

=== ":simple-cmake: CMake"

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

=== ":material-hammer-wrench: Meson"

    ```bash
    cd pineapfel
    meson setup builddir
    ninja -C builddir
    ```

    To install system-wide:

    ```bash
    meson install -C builddir
    ```

!!! note
    If `APFEL++` is installed in a non-standard location, you may need to set `LD_LIBRARY_PATH`
    at runtime:

    ```bash
    export LD_LIBRARY_PATH=/path/to/lib:$LD_LIBRARY_PATH
    ```

---

## Python API

The Python bindings require [pybind11](https://github.com/pybind/pybind11) and
[PineAPPL](https://github.com/NNPDF/pineappl) in addition to the C++ dependencies above.
They are built automatically by Meson when Python 3 and `pybind11` are found.

### Prerequisites

```bash
pip install pybind11 pineappl numpy
```

### Building

=== ":material-hammer-wrench: Meson"

    The Python extension module is compiled alongside the C++ library.  No extra
    configuration flags are needed — Meson detects Python and pybind11 automatically:

    ```bash
    meson setup builddir
    ninja -C builddir
    ```

    A file named `_pineapfel.cpython-<version>-<arch>.so` appears in `builddir/`.

### Installing

=== ":simple-python: System-wide"

    ```bash
    # `pineapfel` is placed in your Python's site-packages alongside the C++ library.
    meson install -C builddir
    ```

=== ":simple-python: In-tree (development)"

    Symlink the built extension into the source package so `import pineapfel` works
    without installation:

    ```bash
    ln -sf "$(pwd)/builddir/_pineapfel.cpython-$(python3 -c 'import sys; v=sys.version_info; print(f"{v.major}{v.minor}")')-"*.so pineapfel_pyapi/pineapfel/

    # Prepend the source tree to PYTHONPATH:
    export PYTHONPATH="$(pwd)/pineapfel_pyapi:$PYTHONPATH"
    python3 -c "import pineapfel; print('OK')"
    ```

!!! note
    See [Python API](python-api.md) for a full usage guide and API reference.
