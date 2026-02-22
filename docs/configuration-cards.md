## Writing configuration cards

PineAPFEL uses YAML configuration files to specify all physics and numerical parameters.
The **theory card** and **operator card** are required for both grid creation and evolution.
An additional **grid card** is needed when creating grids with `build_grid()` or
`pineapfel-build` (see [Grid creation](grid-creation.md) for the full grid card reference).

### Theory card

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

# Physical heavy-quark masses for massive coefficient functions, all 6
# flavors required (up, down, strange, charm, bottom, top) in GeV.
# Used when MassScheme: FFN or MassScheme: FONLL is set in the grid card.
# Defaults to QuarkThresholds padded to 6 entries (top = 172.0 GeV) if absent.
HeavyQuarkMasses: [0.0, 0.0, 0.0, 1.4, 4.75, 172.0]

# CKM squared matrix elements |Vij|^2, 9 entries in row-major order:
# [Vud^2, Vus^2, Vub^2, Vcd^2, Vcs^2, Vcb^2, Vtd^2, Vts^2, Vtb^2]
# Only used for CC grids. Optional: defaults to standard PDG values if absent.
CKM: [0.94922, 0.05077, 0.00001, 0.05073, 0.94760, 0.00168, 0.00007, 0.00162, 0.99831]

# Output flavor PIDs written to the FK table.
# These define the pids_out dimension of the evolution operator.
# Standard PDG codes: 1=d, 2=u, 3=s, 4=c, 5=b, -1=dbar, ..., 21=gluon, 22=photon
Flavors: [-5, -4, -3, -2, -1, 1, 2, 3, 4, 5, 21]

# --- Massive coefficient function tabulation parameters ---
# These control the xi = Q^2/m^2 interpolation grid built by APFEL++ when
# MassScheme: FFN or MassScheme: FONLL is set in the grid card.
# All fields are optional; the values shown are the defaults.

# Number of xi interpolation points.
MassNxi: 150

# xi range [MassXiMin, MassXiMax] over which massive CFs are tabulated.
MassXiMin: 0.05
MassXiMax: 10000.0

# Polynomial interpolation degree in xi.
MassIntDeg: 3

# Lambda parameter for the xi-grid density (controls spacing near xi = 0).
MassLambda: 0.0005

# Interpolation mode passed to TabulateObject (0 = default).
MassIMod: 0

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

### Operator card

```yaml
# x-space interpolation grid.
# Each sub-grid is defined by:
#   n_knots:     number of interpolation knots
#   x_min:       lower edge of the sub-grid in Bjorken x
#   poly_degree: polynomial interpolation degree
#
# The sub-grids are joined into a single APFEL++ Grid object.
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
# Use [1.0, 1.0, 1.0] for central predictions.
xi: [1.0, 1.0, 1.0]
```

