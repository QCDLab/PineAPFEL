## Creating grids

PineAPFEL can create PineAPPL interpolation grids pre-filled with analytically computed
coefficient functions from APFEL++. This is useful for structure function computations
where the hard-scattering kernels are known analytically and only the non-perturbative
input (PDFs or Fragmentation Functions) needs to be provided at convolution time.

The grid creation workflow requires three YAML configuration files:

1. A **grid card** that defines the process, observable, binning, and perturbative orders
2. A **theory card** that specifies the QCD parameters (coupling, thresholds, perturbative order)
3. An **operator card** that defines the x-space interpolation grid and tabulation parameters

### Supported processes and observables

The `build_grid()` function currently supports the following combinations:

| Process | Observable | Current | APFEL++ initializer (ZM) | Mass schemes |
|---------|-----------|---------|--------------------------|-------------|
| DIS | \(F_2\) | NC | `InitializeF2NCObjectsZM` | ZM, FFN, FONLL |
| DIS | \(F_L\) | NC | `InitializeFLNCObjectsZM` | ZM, FFN, FONLL |
| DIS | \(F_3\) | NC | `InitializeF3NCObjectsZM` | ZM only |
| DIS | \(F_2\) | CC\(+\) | `InitializeF2CCPlusObjectsZM` | ZM only |
| DIS | \(F_L\) | CC\(+\) | `InitializeFLCCPlusObjectsZM` | ZM only |
| DIS | \(F_3\) | CC\(+\) | `InitializeF3CCPlusObjectsZM` | ZM only |
| DIS | \(F_2\) | CC\(-\) | `InitializeF2CCMinusObjectsZM` | ZM only |
| DIS | \(F_L\) | CC\(-\) | `InitializeFLCCMinusObjectsZM` | ZM only |
| DIS | \(F_3\) | CC\(-\) | `InitializeF3CCMinusObjectsZM` | ZM only |
| SIA | \(F_2\) | NC | `InitializeF2NCObjectsZMT` | ZM only |
| SIA | \(F_L\) | NC | `InitializeFLNCObjectsZMT` | ZM only |
| SIA | \(F_3\) | NC | `InitializeF3NCObjectsZMT` | ZM only |
| SIDIS | \(F_2\) | NC | `InitializeSIDIS` | ZM only |
| SIDIS | \(F_L\) | NC | `InitializeSIDIS` | ZM only |
| DIS (polarized) | \(g_1\) | NC | `Initializeg1NCObjectsZM` | ZM only |
| DIS (polarized) | \(g_L\) | NC | `InitializegLNCObjectsZM` | ZM only |
| DIS (polarized) | \(g_4\) | NC | `Initializeg4NCObjectsZM` | ZM only |
| SIDIS (polarized) | \(G_1\) | NC | `InitializeSIDISpol` | ZM only |

Polarized grids are selected by setting `Polarized: true` in the grid card. The
`Observable` field retains its unpolarized name (`F2`, `FL`, `F3`) and is interpreted
as the corresponding polarized observable (\(g_1\), \(g_L\), \(g_4\), or \(G_1\) for
SIDIS) when `Polarized: true`.

The mass scheme is selected with the `MassScheme` field in the grid card (see
[Mass scheme](#mass-scheme) below). Combinations that only support ZM will emit a
warning and fall back to ZM when a non-ZM scheme is requested.

!!! info "What is not yet supported"
    - **SIA + CC** is not supported (APFEL++ only provides CC initializers for DIS).
    - **SIDIS + CC** is not supported.
    - **SIDIS + F3** is not supported (APFEL++ does not provide SIDIS F3 coefficient functions).
    - **Polarized + CC** is not supported.
    - **Polarized SIA** is not supported (APFEL++ does not provide time-like polarized SF initializers).
    - **Polarized SIDIS + FL** is not supported (only \(G_1\) is available).
    - **QCD+QED** coefficient functions are not implemented in grid filling (QED corrections
      are only available in the evolution step).

### Convolution types

Each process type requires a specific convolution type, which determines the type of
non-perturbative input the grid will be convoluted with:

| Process | `ConvolutionTypes` | Description |
|---------|-------------------|-------------|
| DIS | `[(UN)POL_PDF]` | (Un)polarised parton distribution functions |
| SIA | `[(UN)POL_FF]` | (Un)polarised fragmentation functions |
| SIDIS | `[(UN)POL_PDF, (UN)POL_FF]` | PDF for the initial state, FF for the final state |

### Mass scheme

The heavy-quark mass scheme is set via the `MassScheme` field in the grid card:

| Value | Description |
|-------|-------------|
| `ZM` | **Zero-mass Variable Flavor Number** `(default)`. Quarks are massless above their thresholds. Applies to all processes and observables. |
| `FFN` | **Fixed-Flavor Number**. Uses `InitializeF2/FLNCObjectsMassive` from APFEL++. Massive charm and bottom coefficient functions depending on \(\xi = Q^2/m^2\) are included. DIS NC \(F_2\)/\(F_L\) only. |
| `FONLL` | **FONLL** scheme. Implemented as \(F_\mathrm{FONLL} = F_\mathrm{ZM} + F_\mathrm{FFN}\), following APFEL++ conventions. Available for DIS NC \(F_2\)/\(F_L\) only. |
| `MassiveZero` | Massless limit of FFN (`InitializeF2/FLNCObjectsMassiveZero`). In APFEL++ the total channel is set to zero for this initializer, so the resulting grid carries zero coefficient functions. Included for completeness. |


!!! warning
    When a non-ZM scheme is requested for an unsupported combination (CC current, \(F_3\),
    polarized, SIA, SIDIS), PineAPFEL prints a warning and falls back to ZM automatically.

The heavy-quark masses used by FFN and FONLL are read from the `HeavyQuarkMasses` field
in the theory card (6 entries, all flavors including top). The ξ-grid tabulation is
controlled by the `MassNxi`, `MassXiMin`, `MassXiMax`, `MassIntDeg`, `MassLambda`, and
`MassIMod` theory-card fields. See [Configuration cards](configuration-cards.md) for
the full list.

### Perturbative orders

#### DIS and SIA (unpolarized)

The coefficient functions are available at the following perturbative orders:

| `alpha_s` power | Label | \(F_2\)/\(F_L\) content | \(F_3\) content |
|:-:|:--:|---|---|
| 0 | LO | \(\delta(1-x)\) for quarks, 0 for gluon | \(\delta(1-x)\) for quarks, 0 for gluon |
| 1 | NLO | \(C_{2,\mathrm{NS}}^{(1)}\), \(C_{2,g}^{(1)}\) | \(C_{3,\mathrm{NS}}^{(1)}\), no gluon |
| 2 | NNLO | \(C_{2,\mathrm{NS}}^{(2)}\), \(C_{2,\mathrm{PS}}^{(2)}\), \(C_{2,g}^{(2)}\) | \(C_{3,\mathrm{NS}}^{(2)}\), no gluon |

#### DIS (polarized)

For polarized DIS (`ConvolutionTypes: [POL_PDF]`), the structure functions \(g_1\), \(g_L\), and \(g_4\)
are selected instead. The available coefficient functions mirror the unpolarized case:

| `alpha_s` power | Label | \(g_1\) / \(g_L\) content | \(g_4\) content |
|:-:|:--:|---|---|
| 0 | LO | \(\delta(1-x)\) for quarks, 0 for gluon | \(\delta(1-x)\) for quarks, 0 for gluon |
| 1 | NLO | \(G_{1,\mathrm{NS}}^{(1)}\), \(G_{1,g}^{(1)}\) | \(G_{4,\mathrm{NS}}^{(1)}\), no gluon |
| 2 | NNLO | \(G_{1,\mathrm{NS}}^{(2)}\), \(G_{1,\mathrm{PS}}^{(2)}\), \(G_{1,g}^{(2)}\) | \(G_{4,\mathrm{NS}}^{(2)}\), no gluon |

#### SIDIS (unpolarized)

SIDIS coefficient functions depend on two momentum-fraction variables (\(x\) and \(z\)) and are
provided as `DoubleObject<Operator>` instances by APFEL++. The available terms per channel type
at each order are:

| `alpha_s` power | Label | \(qq\) | \(gq\) | \(qg\) |
|:-:|:--:|:--:|:--:|:--:|
| 0 | LO | \(C_{2,qq}^{(0)}\) | — | — |
| 1 | NLO | \(C_{2,qq}^{(1)}\) | \(C_{2,gq}^{(1)}\) | \(C_{2,qg}^{(1)}\) |
| 2 | NNLO | \(C_{2,qq}^{(2)}\) (\(n_f\)-dependent) | — | — |

The channel labels refer to the convolution pair (PDF flavour, FF flavour):
\(qq\) = quark PDF ⊗ quark FF, \(gq\) = quark PDF ⊗ gluon FF,
\(qg\) = gluon PDF ⊗ quark FF. The same coefficient functions apply to \(F_L\)
with \(C_{L,qq/gq/qg}\) replacing \(C_{2,\ldots}\); \(F_L\) has no LO contribution.

#### SIDIS (polarized)

For polarized SIDIS (first entry of `ConvolutionTypes` is `POL_PDF`), only \(G_1\) is
available (set `Observable: F2`). The coefficient function structure mirrors the unpolarized
\(F_2\) case:

| `alpha_s` power | Label | \(qq\) | \(gq\) | \(qg\) |
|:-:|:--:|:--:|:--:|:--:|
| 0 | LO | \(G_{1,qq}^{(0)}\) | — | — |
| 1 | NLO | \(G_{1,qq}^{(1)}\) | \(G_{1,gq}^{(1)}\) | \(G_{1,qg}^{(1)}\) |
| 2 | NNLO | \(G_{1,qq}^{(2)}\) (\(n_f\)-dependent) | — | — |

The orders are specified in the grid card via the `Orders` field. Each order entry is a
5-element array `[alpha_s, alpha, log_xir, log_xif, log_xia]`. For pure QCD coefficient
functions, only the `alpha_s` power matters; the remaining entries should be set to 0.

Each entry stores the coefficient function at that **specific** power of \(\alpha_s\),
not the cumulative sum. For a complete NNLO prediction, all three orders (LO, NLO, NNLO)
must be listed so that the grid contains separate subgrids for each perturbative contribution.

!!! warning
    Orders beyond NNLO (`alpha_s > 2`) are silently skipped during grid filling, even
    though APFEL++ provides N3LO coefficient functions for some observables. This is a
    current limitation that may be lifted in a future release.

### Channel decomposition

Channels are **automatically derived** by `build_grid()` from the process, observable, and
the number of active flavours. You do not need to specify them in the grid card. The number
of active flavours \(n_{f}^{\mathrm{max}}\) is determined from the maximum \(Q^2\) across
all bins using the quark thresholds from the theory card.

The `derive_channels()` function generates channels in the **physical (PDG) basis**.

#### DIS and SIA channels

For **\(F_2\)** and **\(F_L\)** (C-even, Neutral Current):

- One quark channel per active flavour \(q = 1, \ldots, n_{f}^{\mathrm{max}}\):
  `pids: [[q], [-q]]`, `factors: [1.0, 1.0]` (i.e. \(q + \bar{q}\))
- One gluon channel: `pids: [[21]]`, `factors: [1.0]`

For **\(F_3\)** (C-odd, Neutral Current):

- One quark channel per active flavour \(q = 1, \ldots, n_{f}^{\mathrm{max}}\):
  `pids: [[q], [-q]]`, `factors: [1.0, -1.0]` (i.e. \(q - \bar{q}\))
- **No gluon channel** (\(C_\mathrm{G} = 0\) at all perturbative orders)

For **Charged-Current (CC)** processes, the channel structure depends on both the
observable and the CC sign variant. The C-parity of the observable determines the
quark combination:

- **C-even** (factors \([1, 1]\), i.e. \(q + \bar{q}\)):
  \(F_2\)/\(F_L\) CC\(+\), \(F_3\) CC\(-\)
- **C-odd** (factors \([1, -1]\), i.e. \(q - \bar{q}\)):
  \(F_2\)/\(F_L\) CC\(-\), \(F_3\) CC\(+\)

A gluon channel is present only for \(F_2/F_L\) with NC or CC\(+\). For CC, the per-quark
weights \(w_q\) are the sum of CKM² elements where quark \(q\) participates
(filtered by active partner flavours), replacing the electroweak charges used in NC.

For example, with 5 active flavours and observable F2 NC, the auto-derived channels are:

| Channel | PIDs | Factors |
|---------|------|---------|
| \(d + \bar{d}\) | `[[1], [-1]]` | `[1.0, 1.0]` |
| \(u + \bar{u}\) | `[[2], [-2]]` | `[1.0, 1.0]` |
| \(s + \bar{s}\) | `[[3], [-3]]` | `[1.0, 1.0]` |
| \(c + \bar{c}\) | `[[4], [-4]]` | `[1.0, 1.0]` |
| \(b + \bar{b}\) | `[[5], [-5]]` | `[1.0, 1.0]` |
| \(g\) | `[[21]]` | `[1.0]` |

The per-channel coefficient functions are constructed from the APFEL++ operators
\(C_\mathrm{NS}\), \(C_\mathrm{S}\), and \(C_\mathrm{G}\) using the general formula:

| Channel | Coefficient function |
|---------|---------------------|
| Quark \(q\) | \(\mathcal{C}_q = w_q \, C_\mathrm{NS} + \frac{\Sigma_w}{6}\,(C_\mathrm{S} - C_\mathrm{NS})\) |
| Gluon | \(\mathcal{C}_g = \Sigma_w \, C_\mathrm{G}\) |

where \(w_q\) is the per-quark weight (electroweak charge \(e_q^2\) for NC, or CKM
weight for CC), \(\Sigma_w = \sum_{i=1}^{n_f^\mathrm{light}} w_i\) sums over light quarks
only, and the factor of 6 matches the internal normalisation convention used in APFEL++'s
`DISNCBasis`/`DISCCBasis`. APFEL++ sets \(C_\mathrm{S} = C_\mathrm{NS}\) and/or
\(C_\mathrm{G} = 0\) where the physics requires it, so the same formula works for all
observables and currents.

For the **FFN and FONLL** schemes, all \(n_f^\mathrm{max}\) quark channels contribute to
the pure-singlet (PS) term \(\frac{\Sigma_w}{6}(C_\mathrm{S}-C_\mathrm{NS})\): light
quarks \((q \leq n_f^\mathrm{light})\) also receive the non-singlet (NS) weight \(w_q\),
while heavy quarks \((q > n_f^\mathrm{light})\) have \(w_q = 0\) for the NS term. This
ensures that the summed PS contribution is proportional to the full SIGMA distribution
\(\sum_{i=1}^{n_f^\mathrm{max}} 2f_i\) as required by APFEL++. Additionally, heavy-quark
massive gluon coefficients \(C_\mathrm{G}^{(\mathrm{h})}\) are accumulated into the gluon
channel weighted by \(e_{q_h}^2\).

#### SIDIS channels

SIDIS grids carry **two convolutions** (PDF ⊗ FF), so each channel entry specifies a pair of
PIDs rather than a single one. Three channel types are generated per active quark
\(q = 1, \ldots, n_{f}^{\mathrm{max}}\):

| Channel type | PIDs | Factors | Description |
|---|---|:---:|---|
| \(qq\) | `[[q, q], [-q, -q]]` | `[1.0, 1.0]` | Quark PDF ⊗ quark FF (and anti-quark) |
| \(gq\) | `[[q, 21], [-q, 21]]` | `[1.0, 1.0]` | Quark PDF ⊗ gluon FF |
| \(qg\) | `[[21, q], [21, -q]]` | `[1.0, 1.0]` | Gluon PDF ⊗ quark FF |

The channel ordering in the grid is `qq`, `gq`, `qg` for quark 1, then `qq`, `gq`, `qg`
for quark 2, and so on. With 5 active flavours this gives 15 channels in total.

The electroweak weight \(e_q^2\) is applied per-quark directly to the subgrid values
during filling; the channel PIDs themselves are charge-neutral in the grid card.

!!! note
    The `Channels` field in the grid card is still accepted for backward compatibility,
    but it is **always overridden** by the auto-derived channels in `build_grid()`.
    It is recommended to omit `Channels` from the grid card entirely.

---

### The grid card

The grid card defines the structure of the PineAPPL grid. Below is a complete reference
with all fields.

```yaml
# The scattering process.
# Supported values: DIS, SIA, SIDIS
Process: DIS

# The structure function observable.
# Supported values: F2, FL, F3
# Optional, defaults to F2.
Observable: F2

# The electroweak current type.
# Supported values: NC, CC
# Optional, defaults to NC.
Current: NC

# CC sign variant (only used when Current: CC).
# Supported values: Plus, Minus
# Plus  = (F(nu) + F(nubar)) / 2
# Minus = (F(nu) - F(nubar)) / 2
# Optional, defaults to Plus.
# CCSign: Plus

# Heavy-quark mass treatment scheme.
# Supported values: ZM, FFN, FONLL, MassiveZero
#   ZM          — zero-mass VFN scheme (default for all processes)
#   FFN         — fixed-flavor-number scheme with massive CFs (DIS NC F2/FL only)
#   FONLL       — FONLL = F_ZM + F_FFN (DIS NC F2/FL only)
#   MassiveZero — massless limit of FFN; outputs zero CFs (APFEL++ convention)
# For unsupported combinations the scheme falls back to ZM with a warning.
# Optional, defaults to ZM.
# MassScheme: ZM

# PID basis for the channel definitions.
# Supported values: PDG, EVOL
PidBasis: PDG

# Hadron PIDs involved in the process (one per convolution).
#   DIS:   [2212]        (proton)
#   SIA:   [211]         (pion)
#   SIDIS: [2212, 211]   (proton + pion)
HadronPids: [2212]

# Convolution types (one per convolution).
#   UNPOL_PDF  — unpolarised parton distributions
#   POL_PDF    — longitudinally polarised parton distributions
#   UNPOL_FF   — unpolarised fragmentation functions
#   POL_FF     — polarised fragmentation functions
#
# Polarization is inferred from the first entry (the PDF slot):
#   POL_PDF / POL_FF  → polarized coefficient functions are used
#     DIS:   F2 -> g1,  FL -> gL,  F3 -> g4  (NC only)
#     SIDIS: F2 -> G1                        (NC only)
#   UNPOL_PDF / UNPOL_FF → unpolarized coefficient functions
# Polarized CC and polarized SIA are not supported.
# For SIDIS the PDF and FF types can differ independently
# (e.g. [POL_PDF, UNPOL_FF] for helicity-weighted cross sections).
ConvolutionTypes: [UNPOL_PDF]

# Perturbative orders as 5-element arrays:
# [alpha_s, alpha, log(xi_R), log(xi_F), log(xi_A)]
# For QCD-only coefficient functions, set all but alpha_s to 0.
Orders:
  - [0, 0, 0, 0, 0]   # O(alpha_s^0) = LO
  - [1, 0, 0, 0, 0]   # O(alpha_s^1) = NLO
  - [2, 0, 0, 0, 0]   # O(alpha_s^2) = NNLO

# Partonic channels (optional).
# Channels are automatically derived by build_grid() from the observable
# and the number of active flavours. You do not need to specify this field.
# If present, it will be overridden during grid building.
#
# Channels:
#   - pids: [[2], [-2]]     # u + ubar
#     factors: [1.0, 1.0]
#   - pids: [[1], [-1]]     # d + dbar
#     factors: [1.0, 1.0]
#   - pids: [[21]]          # gluon
#     factors: [1.0]

# Kinematic bins. Each bin is defined by lower and upper edges
# in each dimension:
#   DIS:   [Q^2, x]        (2 dimensions)
#   SIA:   [Q^2, z]        (2 dimensions)
#   SIDIS: [Q^2, x, z]     (3 dimensions)
Bins:
  - lower: [10.0, 0.001]
    upper: [100.0, 0.01]
  - lower: [100.0, 0.01]
    upper: [1000.0, 0.1]

# Bin normalisation factors (one per bin).
Normalizations: [1.0, 1.0]
```

#### DIS example

A DIS \(F_2\) grid up to NNLO with two \((Q^2, x)\) bins.

```yaml
Process: DIS
Observable: F2
Current: NC
PidBasis: PDG
HadronPids: [2212]
ConvolutionTypes: [UNPOL_PDF]

Orders:
  - [0, 0, 0, 0, 0]
  - [1, 0, 0, 0, 0]
  - [2, 0, 0, 0, 0]

Bins:
  - lower: [10.0, 0.001]
    upper: [100.0, 0.01]
  - lower: [100.0, 0.01]
    upper: [1000.0, 0.1]

Normalizations: [1.0, 1.0]
```

#### SIA example

An SIA \(F_2\) grid for pion production up to NNLO. Note that the second kinematic
dimension is the hadron momentum fraction \(z\) instead of Bjorken \(x\),
and the convolution type is `UNPOL_FF`. Channels are auto-derived as for DIS:

```yaml
Process: SIA
Observable: F2
Current: NC
PidBasis: PDG
HadronPids: [211]
ConvolutionTypes: [UNPOL_FF]

Orders:
  - [0, 0, 0, 0, 0]
  - [1, 0, 0, 0, 0]
  - [2, 0, 0, 0, 0]

Bins:
  - lower: [10.0, 0.2]
    upper: [100.0, 0.4]
  - lower: [100.0, 0.4]
    upper: [1000.0, 0.6]

Normalizations: [1.0, 1.0]
```

#### SIDIS example

A SIDIS \(F_2\) grid for proton→pion semi-inclusive production up to NLO. Bins are
three-dimensional \((Q^2, x, z)\) and two convolution types are required (PDF and FF):

```yaml
Process: SIDIS
Observable: F2
Current: NC
PidBasis: PDG
HadronPids: [2212, 211]
ConvolutionTypes: [UNPOL_PDF, UNPOL_FF]

Orders:
  - [0, 0, 0, 0, 0]   # LO
  - [1, 0, 0, 0, 0]   # NLO

Bins:
  - lower: [10.0, 0.001, 0.2]
    upper: [100.0, 0.01, 0.4]
  - lower: [100.0, 0.01, 0.4]
    upper: [1000.0, 0.1, 0.6]

Normalizations: [1.0, 1.0]
```

#### Polarized DIS example

A polarized DIS \(g_1\) grid up to NLO. Setting `ConvolutionTypes: [POL_PDF]` signals
polarized coefficient functions; `Observable: F2` then selects \(g_1\):

```yaml
Process: DIS
Observable: F2
Current: NC
PidBasis: PDG
HadronPids: [2212]
ConvolutionTypes: [POL_PDF]

Orders:
  - [0, 0, 0, 0, 0]
  - [1, 0, 0, 0, 0]

Bins:
  - lower: [10.0, 0.001]
    upper: [100.0, 0.01]
  - lower: [100.0, 0.01]
    upper: [1000.0, 0.1]

Normalizations: [1.0, 1.0]
```

#### Polarized SIDIS example

A polarized SIDIS \(G_1\) grid up to NLO. Only `Observable: F2` is valid here.
Note that the PDF is polarized while the FF remains unpolarized — this distinction
is expressed directly through `ConvolutionTypes` without any separate flag:

```yaml
Process: SIDIS
Observable: F2
Current: NC
PidBasis: PDG
HadronPids: [2212, 211]
ConvolutionTypes: [POL_PDF, UNPOL_FF]

Orders:
  - [0, 0, 0, 0, 0]
  - [1, 0, 0, 0, 0]

Bins:
  - lower: [10.0, 0.001, 0.2]
    upper: [100.0, 0.01, 0.4]
  - lower: [100.0, 0.01, 0.4]
    upper: [1000.0, 0.1, 0.6]

Normalizations: [1.0, 1.0]
```

#### FFN DIS example

A DIS \(F_2\) grid in the fixed-flavor-number scheme. Requires `HeavyQuarkMasses` in the
theory card (6 entries). The massive coefficient functions are tabulated internally on a
\(\xi = Q^2/m^2\) grid; the tabulation parameters can be tuned with `MassNxi` etc. in
the theory card.

```yaml
Process: DIS
Observable: F2
Current: NC
MassScheme: FFN
PidBasis: PDG
HadronPids: [2212]
ConvolutionTypes: [UNPOL_PDF]

Orders:
  - [0, 0, 0, 0, 0]
  - [1, 0, 0, 0, 0]
  - [2, 0, 0, 0, 0]

Bins:
  - lower: [10.0, 0.001]
    upper: [100.0, 0.01]
  - lower: [100.0, 0.01]
    upper: [1000.0, 0.1]

Normalizations: [1.0, 1.0]
```

#### FONLL DIS example

A DIS \(F_2\) FONLL grid. PineAPFEL computes \(F_\mathrm{FONLL} = F_\mathrm{ZM} + F_\mathrm{FFN}\)
following APFEL++ conventions. The same theory card as for FFN is required.

```yaml
Process: DIS
Observable: F2
Current: NC
MassScheme: FONLL
PidBasis: PDG
HadronPids: [2212]
ConvolutionTypes: [UNPOL_PDF]

Orders:
  - [0, 0, 0, 0, 0]
  - [1, 0, 0, 0, 0]
  - [2, 0, 0, 0, 0]

Bins:
  - lower: [10.0, 0.001]
    upper: [100.0, 0.01]
  - lower: [100.0, 0.01]
    upper: [1000.0, 0.1]

Normalizations: [1.0, 1.0]
```

#### CC DIS example

A DIS \(F_2\) charged-current (CC) grid with the Plus variant \((F(\nu) + F(\bar\nu))/2\).
The `CCSign` field selects between Plus and Minus. The CKM matrix elements are specified
in the theory card (see below):

```yaml
Process: DIS
Observable: F2
Current: CC
CCSign: Plus
PidBasis: PDG
HadronPids: [2212]
ConvolutionTypes: [UNPOL_PDF]

Orders:
  - [0, 0, 0, 0, 0]
  - [1, 0, 0, 0, 0]
  - [2, 0, 0, 0, 0]

Bins:
  - lower: [10.0, 0.001]
    upper: [100.0, 0.01]
  - lower: [100.0, 0.01]
    upper: [1000.0, 0.1]

Normalizations: [1.0, 1.0]
```

The theory card should include a `CKM` field with 9 squared CKM matrix elements
\(|V_{ij}|^2\) in row-major order: \([V_{ud}^2, V_{us}^2, V_{ub}^2, V_{cd}^2,
V_{cs}^2, V_{cb}^2, V_{td}^2, V_{ts}^2, V_{tb}^2]\). If absent, standard PDG
values are used.

---

### Using the CLI

The `pineapfel-build` executable creates a filled PineAPPL grid from three YAML cards:

```text
pineapfel-build <grid.yaml> <theory.yaml> <operator.yaml> [-o output.pineappl.lz4]
```

If `-o` is not specified, the output filename is derived from the grid card by replacing
`.yaml` with `.pineappl.lz4`.

```bash
# Build a DIS F2 grid
pineapfel-build runcards/grid_dis.yaml runcards/theory.yaml runcards/operator.yaml

# Build an SIA grid with a custom output name
pineapfel-build runcards/grid_sia.yaml runcards/theory.yaml runcards/operator.yaml \
    -o sia_f2.pineappl.lz4
```

### Using the library

The same functionality is available programmatically through the `build_grid()` function:

```cpp
#include <pineapfel.h>
#include <pineappl_capi.h>
#include <iostream>

int main() {
    // 1. Load all three cards
    auto grid_def = pineapfel::load_grid_def("runcards/grid_dis.yaml");
    auto theory   = pineapfel::load_theory_card("runcards/theory.yaml");
    auto op_card  = pineapfel::load_operator_card("runcards/operator.yaml");

    // 2. Build and fill the grid
    auto* grid = pineapfel::build_grid(grid_def, theory, op_card);

    // 3. Write the grid
    pineappl_grid_write(grid, "dis_f2.pineappl.lz4");

    // 4. Cleanup
    pineappl_grid_delete(grid);
    return 0;
}
```

The returned grid can also be passed directly to `pineapfel::evolve()` to produce an
FK table in the same program:

```cpp
auto* grid    = pineapfel::build_grid(grid_def, theory, op_card);
auto* fktable = pineapfel::evolve(grid, theory, op_card);
pineappl_grid_write(fktable, "dis_f2.fk.pineappl.lz4");
pineappl_grid_delete(fktable);
pineappl_grid_delete(grid);
```

### Grid structure internals

Understanding the internal layout of the generated grid is useful for debugging and
for writing custom grid-filling code.

#### Node selection

The grid nodes are defined automatically:

- **\(x/z\) nodes**: Taken from the APFEL++ joint interpolation grid, which is built from
  the `xgrid` definition in the operator card. This ensures consistency between the
  coefficient function grid and any subsequent evolution step.

- **\(Q^2\) nodes**: Derived from the bin edges in the grid card. The bin boundaries in the
  first dimension (\(Q^2\)) are collected, and geometrically spaced intermediate points are
  added within each bin. Additionally, quark threshold values (\(m_q^2\)) that fall within
  the overall \(Q^2\) range are included to capture flavour-threshold effects.

#### Subgrid layout

##### DIS and SIA

Each subgrid (one per combination of bin, perturbative order, and channel) is a
two-dimensional array of shape `[n_Q2, n_x]`, stored in row-major order. The
`node_values` vector concatenates the \(Q^2\) nodes followed by the \(x\)/\(z\) nodes:

```
node_values = [Q^2_0, ..., Q^2_{nq-1}, x_0, ..., x_{nx-1}]
```

For each \(Q^2\) node, the coefficient function operator is evaluated at the bin's
\(x\)/\(z\) centre (geometric mean of the bin edges) to produce a distribution on the
APFEL++ joint grid, which populates one row of the subgrid.

##### SIDIS

SIDIS subgrids are three-dimensional arrays of shape `[n_Q2, n_x, n_z]`, stored in
row-major order. The same APFEL++ joint grid is used for both the \(x\) (PDF) and
\(z\) (FF) dimensions. The `node_values` vector concatenates three segments:

```
node_values = [Q^2_0, ..., Q^2_{nq-1}, x_0, ..., x_{nx-1}, z_0, ..., z_{nz-1}]
```

For each \(Q^2\) node and each `DoubleObject` term in the coefficient function,
the \(x\)-distribution is evaluated at the bin's \(x\)-centre and the \(z\)-distribution
at the bin's \(z\)-centre (both geometric means). The subgrid entry at `[iq, ix, iz]`
is the outer product of these two distributions, weighted by \(e_q^2\) and the term
coefficient:

```
subgrid[iq, ix, iz] = e_q^2 * sum_terms { c_term * K_x(x[ix]; x_centre) * K_z(z[iz]; z_centre) }
```

where \(K_x\) and \(K_z\) are the APFEL++ distribution kernels in the \(x\) and \(z\)
directions respectively.

### Programmatic grid definition

In addition to loading from YAML, you can construct a `GridDef` programmatically.
Note that `channels` can be left empty — `build_grid()` will auto-derive them:

```cpp
pineapfel::GridDef def;
def.process    = pineapfel::ProcessType::DIS;
def.observable = pineapfel::Observable::F2;
def.current    = pineapfel::Current::NC;
def.pid_basis  = PINEAPPL_PID_BASIS_PDG;
def.hadron_pids        = {2212};
def.convolution_types  = {PINEAPPL_CONV_TYPE_UNPOL_PDF};

def.orders = {{0, 0, 0, 0, 0}, {1, 0, 0, 0, 0}, {2, 0, 0, 0, 0}};  // LO + NLO + NNLO
// channels are auto-derived by build_grid() — no need to set them
def.bins = {
    {{10.0, 0.001}, {100.0, 0.01}},
    {{100.0, 0.01}, {1000.0, 0.1}},
};
def.normalizations = {1.0, 1.0};

auto* grid = pineapfel::build_grid(def, theory, op_card);
```

You can also call `derive_channels()` directly if you need the channels before
calling `build_grid()`:

```cpp
// Derive channels for DIS F2 NC with 5 active flavours
auto channels = pineapfel::derive_channels(
    pineapfel::ProcessType::DIS,
    pineapfel::Observable::F2,
    pineapfel::Current::NC,
    pineapfel::CCSign::Plus,
    5);
// Returns 6 channels: d+dbar, u+ubar, s+sbar, c+cbar, b+bbar, gluon

// Derive channels for SIDIS F2 with 5 active flavours
auto sidis_channels = pineapfel::derive_channels(
    pineapfel::ProcessType::SIDIS,
    pineapfel::Observable::F2,
    pineapfel::Current::NC,
    pineapfel::CCSign::Plus,
    5);
// Returns 15 channels: (qq, gq, qg) x 5 quarks
```
