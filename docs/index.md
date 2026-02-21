---
title: ""
hide:
  - title
---

<div style="display: flex; justify-content: center;">
  <img src="./assets/pineapfel.svg" width="650" />
</div>

PineAPFEL is a C++ library and CLI tool for creating and evolving [PineAPPL](https://github.com/NNPDF/pineappl)
interpolation grids using [APFEL++](https://github.com/vbertone/apfelxx). It can fill grids with analytically
computed coefficient functions for unpolarized and longitudinally polarized DIS, SIA, and SIDIS structure
functions, and evolve them into FK tables via DGLAP evolution.

All physics parameters (coupling constants, thresholds, perturbative order, etc.) are specified through YAML
configuration files — a **theory card**, an **operator card**, and optionally a **grid card** — rather than
being hardcoded.

## Get started

- [Building](building.md)
- [Quick start](quick-start.md)

## Usage

- [Configuration cards](configuration-cards.md)
- [Grid creation](grid-creation.md)
- [CLI](cli.md)
- [Library](library.md)
- [QCD vs QCD⊗QED](evolution-modes.md)

