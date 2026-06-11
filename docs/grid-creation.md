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
| SIDIS | \(F_T\) | NC | `InitializeSidisObjects` | ZM only |
| SIDIS | \(F_L\) | NC | `InitializeSidisObjects` | ZM only |
| DIS (polarized) | \(g_1\) | NC | `Initializeg1NCObjectsZM` | ZM only |
| DIS (polarized) | \(g_L\) | NC | `InitializegLNCObjectsZM` | ZM only |
| DIS (polarized) | \(g_4\) | NC | `Initializeg4NCObjectsZM` | ZM only |
| SIDIS (polarized) | \(G_1\) | NC | `InitializeSidisObjects` | ZM only |

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
provided as `DoubleOperator` instances by APFEL++ via `InitializeSidisObjects()`. This is the
**exact NNLO** implementation (arXiv:2401.16281), covering all nine partonic channel types.
The available channels per perturbative order are:

| `alpha_s` power | Label | \(qq\) | \(gq\) | \(qg\) | \(gg\) | \(ps\) | \(q\bar{q}\) | \(qpq\) |
|:-:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|
| 0 | LO | \(C_{T,qq}^{(0)}\) | — | — | — | — | — | — |
| 1 | NLO | \(C_{T,qq}^{(1)}\) | \(C_{T,gq}^{(1)}\) | \(C_{T,qg}^{(1)}\) | — | — | — | — |
| 2 | NNLO | \(C_{T,\mathrm{NS}}^{(2)}\) (\(n_f\)-dep.) | \(C_{T,gq}^{(2)}\) | \(C_{T,qg}^{(2)}\) | \(C_{T,gg}^{(2)}\) | \(C_{T,\mathrm{PS}}^{(2)}\) | \(C_{T,q\bar{q}}^{(2)}\) | \(C_{T,qpq_{1,2,3}}^{(2)}\) |

The same structure applies to \(F_L\) (replacing \(C_T\) with \(C_L\)); \(F_L\) has no LO
contribution.

The nine NNLO channel labels refer to the convolution pair (PDF flavour, FF flavour):

| Channel | PDF side | FF side | Notes |
|---------|----------|---------|-------|
| \(qq\) (NS) | quark \(q\) | quark \(q\) | \(n_f\)-dependent non-singlet |
| \(gq\) | quark \(q\) | gluon | |
| \(qg\) | gluon | quark \(q\) | |
| \(gg\) | gluon | gluon | |
| \(ps\) | quark \(q\) | quark \(q\) | pure-singlet, weight \(\sum_i e_i^2\) |
| \(q\bar{q}\) | quark \(q\) | antiquark \(\bar{q}\) | |
| \(qpq_1\) | quark \(j\) | quark \(k \neq j\) | sum over target FF |
| \(qpq_2\) | quark \(k \neq i\) | quark \(i\) | sum over source PDF |
| \(qpq_3\) | quark \(a\) | quark \(b \neq a\) | charge product \(e_a e_b\) weight |

#### SIDIS (polarized)

For polarized SIDIS (first entry of `ConvolutionTypes` is `POL_PDF`), only \(G_1\) is
available (set `Observable: F2`). The exact NNLO implementation (arXiv:2404.08597) covers
the same nine channel types. Note that at NNLO the \(gq\) and \(qg\) channels become
\(n_f\)-dependent for the polarized case:

| `alpha_s` power | Label | \(qq\) | \(gq\) | \(qg\) | \(gg\) | \(ps\) | \(q\bar{q}\) | \(qpq\) |
|:-:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|
| 0 | LO | \(G_{1,qq}^{(0)}\) | — | — | — | — | — | — |
| 1 | NLO | \(G_{1,qq}^{(1)}\) | \(G_{1,gq}^{(1)}\) | \(G_{1,qg}^{(1)}\) | — | — | — | — |
| 2 | NNLO | \(G_{1,\mathrm{NS}}^{(2)}\) (\(n_f\)-dep.) | \(G_{1,gq}^{(2)}\) (\(n_f\)-dep.) | \(G_{1,qg}^{(2)}\) (\(n_f\)-dep.) | \(G_{1,gg}^{(2)}\) | \(G_{1,\mathrm{PS}}^{(2)}\) | \(G_{1,q\bar{q}}^{(2)}\) | \(G_{1,qpq_{1,2,3}}^{(2)}\) |

The orders are specified in the grid card via the `Orders` field. Each order entry is a
5-element array `[alpha_s, alpha, log_xir, log_xif, log_xia]`. For central-scale QCD
coefficient functions set all but `alpha_s` to 0. To include renormalization-scale
logarithms (DIS/SIA only), set `log_xir` to the desired power; see
[Renormalization scale variation](#renormalization-scale-variation-dissia) below.

Each entry stores the coefficient function at that **specific** power of \(\alpha_s\)
(and \(\ln\xi_R^2\)), not the cumulative sum. For a complete NNLO prediction, all three
orders (LO, NLO, NNLO) must be listed so that the grid contains separate subgrids for
each perturbative contribution.

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
PIDs rather than a single one. The complete set of channel types depends on the perturbative
order requested.

**LO and NLO** (3 types × \(n_f\) quarks):

| Channel type | PIDs | Factors | Description |
|---|---|:---:|---|
| \(qq\) | `[[q, q], [-q, -q]]` | `[1.0, 1.0]` | Quark PDF ⊗ quark FF (and anti-quark) |
| \(gq\) | `[[q, 21], [-q, 21]]` | `[1.0, 1.0]` | Quark PDF ⊗ gluon FF |
| \(qg\) | `[[21, q], [21, -q]]` | `[1.0, 1.0]` | Gluon PDF ⊗ quark FF |

**Additional NNLO-only channels** (shared across all quarks):

| Channel type | PIDs | Factors | Description |
|---|---|:---:|---|
| \(gg\) | `[[21, 21]]` | `[1.0]` | Gluon PDF ⊗ gluon FF; weight \(\sum_i e_i^2\) |
| \(ps\) (per quark \(q\)) | `[[q, q], [-q, -q]]` | `[1.0, 1.0]` | Pure-singlet; weight \(\sum_i e_i^2\) |
| \(q\bar{q}\) (per quark \(q\)) | `[[q, -q], [-q, q]]` | `[1.0, 1.0]` | Quark PDF ⊗ antiquark FF |
| \(qpq_1\) (per source \(j\)) | `[[j,k],[-j,k],[j,-k],[-j,-k]]` | `[1.0,...]` | Sum over \(k \neq j\); weight \(e_j^2\) |
| \(qpq_2\) (per target \(i\)) | `[[k,i],[k,-i],[-k,i],[-k,-i]]` | `[1.0,...]` | Sum over \(k \neq i\); weight \(e_i^2\) |
| \(qpq_3\) (per pair \(a \neq b\)) | `[[a,b],[-a,b],[a,-b],[-a,-b]]` | `[1.,-1.,-1.,1.]` | Weight \(e_a \cdot e_b\) (may be negative) |

Channels are **always auto-derived** from \(n_f^{\mathrm{max}}\) regardless of the orders
requested — the fill loop simply returns zero for channel–order combinations where no
coefficient function exists (e.g. \(gg\) at LO or NLO).

The channel ordering in the grid is: `qq`, `gq`, `qg` per quark (1…\(n_f\)), then `gg`,
`ps` per quark, `qbq` per quark, `qpq1` per quark, `qpq2` per quark, and finally
`qpq3` for each ordered pair \((a, b)\) with \(a \neq b\).
With 3 active flavours this gives \(3 \times 3 + 1 + 3 + 3 + 3 + 3 + 6 = 22\) channels in total.

The electroweak weight (\(e_q^2\), \(\sum_i e_i^2\), or \(e_a e_b\)) is applied per-channel
directly to the subgrid values during filling; the channel PIDs themselves are charge-neutral
in the grid card.

!!! note
    The `Channels` field in the grid card is still accepted for backward compatibility,
    but it is **always overridden** by the auto-derived channels in `build_grid()`.
    It is recommended to omit `Channels` from the grid card entirely.

---

### Operator card: SIDIS-only settings (`sidis_*` keys) {#sidis-operator-card}

Two optional fields in the **operator card** control SIDIS-specific behaviour.
They are ignored for DIS and SIA.

| Field | Default | Purpose |
|-------|---------|---------|
| `sidis_mode` | `legacy` | `legacy` = Gauss–Legendre \(z\) quadrature inside each \(z\) bin; `bsf_exact` = APFEL++-style interpolation weights (`IntInterpolant`) over the \(z\) range. |
| `sidis_int_eps` | `1e-3` | Relative tolerance for `InitializeSidisObjects` (SIDIS NNLO coefficient-function setup). Increase to `1e-1` for faster but coarser initialisation in tests. |

For APFEL++-aligned (BSF-style) SIDIS grids, set `sidis_mode: bsf_exact` explicitly
in the operator YAML:

```yaml
sidis_mode: bsf_exact
sidis_int_eps: 1e-3   # optional; 1e-1 is enough for quick tests
```

See also [Configuration cards — operator card](configuration-cards.md#operator-card-sidis-fields).

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

# Kinematic points. Each entry is a coordinate tuple:
#   DIS/SIA:  [Q^2, x]              — fully pointwise
#   SIDIS:    [Q^2, x, z_lo, z_hi]  — Q^2 and x are pointwise; z is an integration range
Points:
  - [10.0, 0.001]
  - [100.0, 0.01]

# Bin normalisation factors (one per bin). Optional — defaults to 1.0 per bin.
# Normalizations: [1.0, 1.0]

# Legacy alternative to Points (still accepted for backward compatibility).
# Each bin is defined by lower and upper edges in each dimension:
#   DIS/SIA:  [Q^2, x]     (2 dimensions)
#   SIDIS:    [Q^2, x, z]  (3 dimensions)
# The Q^2 node is taken as the geometric centre sqrt(Q^2_lo * Q^2_hi).
# Bins:
#   - lower: [10.0, 0.001]
#     upper: [100.0, 0.01]

# ── FixConvolutions (optional, SIDIS only) ──────────────────────────────────
#
# Fold one or more convolution slots into the grid at build time by evaluating
# a NeoPDF set at each node.  Each entry reduces the number of convolutions in
# the output grid by one; entries are applied in the order listed.
#
# Fields:
#   index      — 0-based index of the convolution slot to fix.  For a SIDIS
#                grid [UNPOL_PDF, UNPOL_FF] use 0 to fix the PDF or 1 to fix
#                the FF.
#   pdf_set    — NeoPDF set name (same identifier used by neopdf_pdf_load).
#   pdf_member — member index within the set (default 0).
#   xi         — factorization-scale variation factor forwarded to PineAPPL
#                (default 1.0 = central scale).
#
# After all fix steps the output grid has (N − k) convolutions, where k is the
# number of entries.  It can then be evolved with pineapfel-evolve like any
# single-convolution grid.
#
# FixConvolutions:
#   - index: 0
#     pdf_set: NNPDF40_nnlo_as_01180
#     pdf_member: 0    # optional
#     xi: 1.0          # optional
```

#### DIS example

A DIS \(F_2\) grid up to NNLO at two \((Q^2, x)\) points.

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

Points:
  - [10.0,   0.001]
  - [1000.0, 0.1]
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

Points:
  - [10.0,   0.2]
  - [1000.0, 0.6]
```

#### SIDIS example

A SIDIS \(F_T\) grid for proton→pion semi-inclusive production up to NLO. Each point
specifies \((Q^2, x)\) exactly and a \(z\) integration range \([z_\mathrm{lo},
z_\mathrm{hi}]\). Two convolution types are required (PDF and FF):

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

# [Q^2, x, z_lo, z_hi] — Q^2 and x are pointwise; z is an integration range.
Points:
  - [10.0,   0.001, 0.2, 0.4]
  - [1000.0, 0.1,   0.4, 0.6]
```

#### SIDIS NNLO example

Adding NNLO extends the channel set from 3 per quark to 9 types (see
[SIDIS channels](#sidis-channels)). Simply add the NNLO order entry — channels
are auto-derived and all 9 NNLO channel types are filled automatically:

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
  - [2, 0, 0, 0, 0]   # NNLO (exact, arXiv:2401.16281)

Points:
  - [10.0,   0.001, 0.2, 0.4]
  - [1000.0, 0.1,   0.4, 0.6]
```

!!! note
    The NNLO SIDIS coefficient functions are computed exactly from the full 2D
    `DoubleExpression` classes in APFEL++ via `InitializeSidisObjects()`. This replaces
    the former approximated NNLO (non-singlet \(qq\) only, based on threshold
    resummation). The exact implementation covers all nine partonic channels for both
    the unpolarised \(F_T\)/\(F_L\) (arXiv:2401.16281) and the polarised \(G_1\)
    (arXiv:2404.08597).

#### Fixing a convolution with a PDF/FF set {#fix-convolutions}

The optional `FixConvolutions` key in the grid card instructs `pineapfel-build` to fold
one or more convolution slots into the grid immediately after filling it with coefficient
functions, and **before** the grid is written to disk.

This is primarily useful for SIDIS, which naturally produces a two-convolution grid
(PDF ⊗ FF). Fixing the PDF slot yields a single-convolution grid that depends only on
the FF; the result can be evolved with `pineapfel-evolve` in the same way as any DIS or
SIA single-convolution grid.

**Grid card with `FixConvolutions`:**

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

Points:
  - [10.0,   0.001, 0.2, 0.4]
  - [1000.0, 0.1,   0.4, 0.6]

# Fix convolution slot 0 (UNPOL_PDF) with a NeoPDF set.
# The written grid has a single convolution (UNPOL_FF only).
FixConvolutions:
  - index: 0
    pdf_set: NNPDF40_nnlo_as_01180
    pdf_member: 0   # optional, default 0
    xi: 1.0         # optional, default 1.0
```

The `FixConvolutions` list is ordered: if two entries are given the second entry's
`index` refers to the convolution numbering **after** the first fix has been applied.

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

Points:
  - [10.0,   0.001]
  - [1000.0, 0.1]
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

Points:
  - [10.0,   0.001, 0.2, 0.4]
  - [1000.0, 0.1,   0.4, 0.6]
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

Points:
  - [10.0,   0.001]
  - [1000.0, 0.1]
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

Points:
  - [10.0,   0.001]
  - [1000.0, 0.1]
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

Points:
  - [10.0,   0.001]
  - [1000.0, 0.1]
```

The theory card should include a `CKM` field with 9 squared CKM matrix elements
\(|V_{ij}|^2\) in row-major order: \([V_{ud}^2, V_{us}^2, V_{ub}^2, V_{cd}^2,
V_{cs}^2, V_{cb}^2, V_{td}^2, V_{ts}^2, V_{tb}^2]\). If absent, standard PDG
values are used.

---

### Renormalization scale variation {#renormalization-scale-variation-dissia}

PineAPFEL supports renormalization-scale variation for **DIS and SIA** grids. When one or
more orders with `log_xir > 0` are requested, `build_grid()` automatically computes the
corresponding coefficient-function contributions and stores them as separate PineAPPL
subgrids. At convolution time, PineAPPL evaluates

$$
F(x, Q^2;\, \xi_R) \;=\;
\sum_{n,\,m} \alpha_s(\mu_R)^n \,
\bigl[\ln\xi_R^2\bigr]^m \,
W_{n,m}
$$

where \(\xi_R = \mu_R / Q\) and \(W_{n,m}\) are the stored subgrid weights for order
`[n, 0, m, 0, 0]`. Setting \(\xi_R = 1\) recovers the central-scale result; varying it
around 1 gives the renormalization-scale uncertainty band.

#### Available renorm-scale orders

The following `log_xir` orders are derived automatically from the central-scale
coefficient functions \(C_0\), \(C_1\), \(C_2\) via the renormalization group:

| Order entry | Label | Stored weight | Derivation |
|-------------|-------|---------------|------------|
| `[1, 0, 1, 0, 0]` | NLO × \(\ln\xi_R^2\) | \(\beta_0(n_f)\,\tfrac{1}{4\pi}\,C_0\) | \(\partial_{\ln\mu_R^2} F\big\lvert_{\mathrm{NLO}}\) |
| `[2, 0, 1, 0, 0]` | NNLO × \(\ln\xi_R^2\) | \(\beta_0(n_f)\,\tfrac{1}{4\pi}\,C_1\) | \(\partial_{\ln\mu_R^2} F\big\lvert_{\mathrm{NNLO}}\) |
| `[2, 0, 2, 0, 0]` | NNLO × \(\ln^2\xi_R^2\) | \(\tfrac{\beta_0(n_f)^2}{2}\,\tfrac{1}{(4\pi)^2}\,C_0\) | \(\tfrac{1}{2}\partial^2_{\ln\mu_R^2} F\big\lvert_{\mathrm{NNLO}}\) |

Here \(\beta_0(n_f) = 11 - \tfrac{2}{3}n_f\) (the one-loop QCD \(\beta\)-function
coefficient) and the factors of \(\tfrac{1}{4\pi}\) absorb the difference between
APFEL++'s \((\alpha_s/4\pi)^n\) convention and PineAPPL's \(\alpha_s^n\) convention,
\(n_f\) is the number of active flavours at the Q² node.

The central-scale orders `[0,0,0,0,0]`, `[1,0,0,0,0]`, `[2,0,0,0,0]` must also be
included in `Orders` for the corresponding \(C_0\), \(C_1\), \(C_2\) subgrids to exist
(scale-log orders are derived from them). If a base order is absent, the derived log
order will be empty.

#### Example: DIS NLO + NNLO with renorm-scale variation

```yaml
Process: DIS
Observable: F2
Current: NC
PidBasis: PDG
HadronPids: [2212]
ConvolutionTypes: [UNPOL_PDF]

Orders:
  - [0, 0, 0, 0, 0]   # LO           (C_0)
  - [1, 0, 0, 0, 0]   # NLO          (C_1)
  - [2, 0, 0, 0, 0]   # NNLO         (C_2)
  - [1, 0, 1, 0, 0]   # NLO  × ln ξ_R²
  - [2, 0, 1, 0, 0]   # NNLO × ln ξ_R²
  - [2, 0, 2, 0, 0]   # NNLO × ln² ξ_R²

Points:
  - [10.0, 0.001]
```

The six-subgrid grid can then be convoluted at any \(\xi_R\) by passing a non-unity
scale factor to `pineappl_grid_convolve_with_one` (or the Python `Grid.convolve`
wrapper).

!!! warning "DIS/SIA only — SIDIS and factorization logs not yet implemented"
    Renormalization-scale logs are currently filled only for **DIS and SIA** processes.
    SIDIS grids ignore `log_xir > 0` entries (no error is raised; the subgrid is simply
    left empty). **Factorization-scale logs** (`log_xif > 0`) require convolving with
    DGLAP splitting functions and are not yet implemented for any process.

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

#### SIDIS with a fixed PDF convolution {#fix-convolutions-cli}

When the grid card contains `FixConvolutions`, `pineapfel-build` automatically evaluates
the specified NeoPDF set and folds it into the grid before writing, producing a
single-convolution output (FF only in this example). The resulting grid can be fed
directly to `pineapfel-evolve`:

```bash
# Build the SIDIS grid and fix the PDF slot in one step
pineapfel-build runcards/grid_sidis_fixed_pdf.yaml runcards/theory.yaml runcards/operator.yaml -o sidis_ff_only.pineappl.lz4

# Evolve the FF-only grid into an FK table (single-convolution, FF evolution)
pineapfel-evolve sidis_ff_only.pineappl.lz4 runcards/theory.yaml runcards/operator.yaml -o sidis_ff_only.fk.pineappl.lz4
```

Without `FixConvolutions` the two-convolution SIDIS grid is written as-is and requires
both a PDF and an FF when convolving:

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

- **\(Q^2\) nodes (DIS and SIA)**: Each bin contributes exactly **one** \(Q^2\) node,
  taken from the `Points` entry directly. All coefficient functions are evaluated at
  that single \(Q^2\) value, and the subgrid has shape `[1, n_x]`.

- **\(Q^2\) nodes (SIDIS)**: Same pointwise strategy — each bin contributes exactly one
  \(Q^2\) node equal to the \(Q^2\) coordinate from the `Points` entry (geometric centre
  \(\sqrt{Q^2_\mathrm{lo} \cdot Q^2_\mathrm{hi}}\) when the old `Bins` format is used).
  The subgrid has shape `[1, n_x, n_z]`.

#### Subgrid layout

##### DIS and SIA

Each subgrid (one per combination of bin, perturbative order, and channel) is a
two-dimensional array of shape `[1, n_x]`, stored in row-major order. The
`node_values` vector is:

```text
node_values = [Q^2, x_0, ..., x_{nx-1}]
```

The coefficient function operator is evaluated at the bin's \(Q^2\) and \(x\) point to
produce a distribution on the APFEL++ joint grid, which fills the single row.

##### SIDIS

SIDIS subgrids are three-dimensional arrays of shape `[1, n_x, n_z]`, stored in
row-major order. The same APFEL++ joint grid is used for both the \(x\) (PDF) and
\(z\) (FF) dimensions. The `node_values` vector is:

```text
node_values = [Q^2, x_0, ..., x_{nx-1}, z_0, ..., z_{nz-1}]
```

Each coefficient function is stored as a `DoubleOperator` — a full 2D kernel on the
\((x, z)\) grid. The 2D operator is evaluated at the bin's \((Q^2, x)\) point and over
the bin's \(z\) range \([z_\mathrm{lo}, z_\mathrm{hi}]\) using
`eval_double_op_column()`. The subgrid entry at `[ix, iz]` accumulates the weighted
kernel:

```text
subgrid[ix, iz] = fill_weight * W(x[ix], z[iz]; x_point, z_centre)
```

where `fill_weight` is the per-channel electroweak factor (\(e_q^2\), \(\sum_i e_i^2\),
or \(e_a e_b\) depending on channel type) and \(W\) is the 2D interpolation kernel
extracted from the `DoubleOperator`. The \(z\) direction accumulates contributions
from APFEL++-style interpolation weights over joint-grid nodes (`sidis_mode: bsf_exact`)
or from internal Gauss–Legendre quadrature on the \(z\) range (`sidis_mode: legacy`).
The same DoubleOperator column machinery is used for all perturbative orders (LO, NLO, NNLO).

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
// Pointwise bins: lower == upper == {Q^2, x}
def.bins = {
    {{10.0,   0.001}, {10.0,   0.001}},
    {{1000.0, 0.1},   {1000.0, 0.1}},
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

// Derive channels for SIDIS F2 with 5 active flavours (includes all NNLO types)
auto sidis_channels = pineapfel::derive_channels(
    pineapfel::ProcessType::SIDIS,
    pineapfel::Observable::F2,
    pineapfel::Current::NC,
    pineapfel::CCSign::Plus,
    5);
// Returns 50 channels:
//   15 LO/NLO  : (qq, gq, qg) x 5 quarks
//    1 NNLO gg
//    5 NNLO ps  (one per quark)
//    5 NNLO qbq (one per quark)
//    5 NNLO qpq1 (one per source quark)
//    5 NNLO qpq2 (one per target quark)
//   20 NNLO qpq3 (one per ordered pair a≠b, 5x4=20)
// Channels that have no coefficient function at a given order contribute zero.
```
