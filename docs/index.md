---
title: ""
---

<div class="home-page-logo" style="display: flex; justify-content: center;">
  <img src="./assets/pineapfel.svg" width="650" />
</div>

`PineAPFEL` is a C++ library and CLI tool for creating and evolving [PineAPPL](https://github.com/NNPDF/pineappl)
interpolation grids using [APFEL++](https://github.com/vbertone/apfelxx). It can
fill grids with analytically computed coefficient functions for unpolarized and
longitudinally polarized DIS, SIA, and SIDIS structure functions, and evolve them
into FK tables via DGLAP evolution.

All physics parameters (coupling constants, thresholds, perturbative order, etc.)
are specified through YAML configuration files through a **theory card**, an
**operator card**, and optionally a **grid card** rather than being hardcoded.

## Getting started

- [Building & Installing](building.md)
- [Writing configuration cards](configuration-cards.md)
- [Producing PineAPPL grids](grid-creation.md)
- [PineAPFEL CLI](cli.md)
- [Using PineAPFEL as a Library](library.md)
