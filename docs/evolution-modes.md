The evolution mode is controlled by the `QED` field in the theory card:

| Mode | `QED` | Basis size | DGLAP objects | Coupling |
|---|---:|---:|---|---|
| Pure QCD | `false` | 13 × 13 | `apfel::DglapObjects` | `apfel::AlphaQCD` |
| QCD+QED | `true` | 20 × 20 | `apfel::DglapObjectsQCDQED` | `apfel::AlphaQCDQED` |

In **QCD mode**, PineAPFEL supports three convolution types (determined automatically from the grid):

- `UNPOL_PDF` — unpolarised PDFs
- `POL_PDF` — longitudinally polarised PDFs
- `UNPOL_FF` — unpolarised fragmentation functions

In **QCD+QED mode**, PineAPFEL currently supports:

- `UNPOL_PDF` — unpolarised PDFs with QED corrections
