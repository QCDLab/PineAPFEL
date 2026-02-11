## Writing configuration cards

All physics and numerical parameters are split between two YAML files: a **theory card** and an **operator card**.

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

