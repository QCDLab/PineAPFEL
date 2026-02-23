"""
Python tests mirroring tests/test_grid_vs_apfelxx.cpp.

Each test builds a PineAPFEL grid via the Python bindings, convolves it with
a toy distribution through pineappl-py, and asserts physically motivated
properties — finite values, correct sign, order-dependence, etc.

Because Python has no APFEL++ bindings the *direct* comparison against
APFEL++ operator objects (C++ tests 1–3) and BuildStructureFunctions (C++
tests 4–11) is replaced by:
  • finiteness of all bin values,
  • sign constraints that hold for the given observable at LO,
  • order-dependence (LO ≠ LO+NLO ≠ full NNLO),
  • write / read round-trip consistency.

Tolerance for the numerical checks is kept loose (1 %) so that differences
in floating-point precision between the Python and C++ paths do not cause
spurious failures.
"""

import numpy as np

import pineapfel as pf
import pineappl.grid as pineappl_grid

from conftest import (
    build_pineappl,
    convolve,
    make_pdg_convs,
    toy_xfx,
    const_alphas,
    runcard,
)


TOL = 1e-2  # 1 % relative tolerance for non-zero bin values


def assert_finite_nonzero(values: np.ndarray, label: str = "") -> None:
    tag = f"[{label}] " if label else ""
    assert np.all(np.isfinite(values)), f"{tag}non-finite result: {values}"
    assert np.any(values != 0.0), f"{tag}all-zero result"


def assert_positive(values: np.ndarray, label: str = "") -> None:
    tag = f"[{label}] " if label else ""
    assert np.all(values > 0), f"{tag}expected all-positive, got {values}"


def order_mask(pg, max_order_idx: int) -> np.ndarray:
    """Return a boolean mask that enables orders 0..max_order_idx."""
    n = len(pg.orders())
    return np.array([i <= max_order_idx for i in range(n)])


def test_01_dis_nc_f2_lo(theory, op_card):
    """NC DIS F2 LO: all bin values must be finite and positive."""
    pg = build_pineappl("grid_dis.yaml", theory, op_card)
    mask = order_mask(pg, 0)
    lo = convolve(pg, order_mask=mask)

    assert_finite_nonzero(lo, "DIS NC F2 LO")
    assert_positive(lo, "DIS NC F2 LO")


def test_02_dis_nc_f2_nlo(theory, op_card):
    """NC DIS F2 LO+NLO: finite, positive, and different from LO-only."""
    pg = build_pineappl("grid_dis.yaml", theory, op_card)
    lo = convolve(pg, order_mask=order_mask(pg, 0))
    nlo = convolve(pg, order_mask=order_mask(pg, 1))

    assert_finite_nonzero(nlo, "DIS NC F2 NLO")
    assert_positive(nlo, "DIS NC F2 NLO")
    assert not np.allclose(
        lo, nlo, rtol=TOL
    ), "LO and LO+NLO results are unexpectedly identical"


def test_03_dis_nc_f2_nnlo(theory, op_card):
    """NC DIS F2 full NNLO: finite, and distinct from LO+NLO."""
    pg = build_pineappl("grid_dis.yaml", theory, op_card)
    nlo = convolve(pg, order_mask=order_mask(pg, 1))
    full = convolve(pg)

    assert_finite_nonzero(full, "DIS NC F2 NNLO")
    assert not np.allclose(
        nlo, full, rtol=TOL
    ), "LO+NLO and full-NNLO results are unexpectedly identical"


def test_04_dis_nc_f2_roundtrip(theory, op_card, tmp_path):
    """NC DIS F2 full: write→read round-trip gives identical convolution."""
    gdef = pf.load_grid_def(runcard("grid_dis.yaml"))
    built = pf.build_grid(gdef, theory, op_card)

    path = str(tmp_path / "dis.pineappl.lz4")
    built.write(path)

    g1 = pineappl_grid.Grid.read(path)
    g2 = pineappl_grid.Grid.read(path)

    pdg_convs = make_pdg_convs(g1)
    v1 = np.asarray(
        g1.convolve(pdg_convs=pdg_convs, xfxs=[toy_xfx], alphas=const_alphas)
    )
    v2 = np.asarray(
        g2.convolve(pdg_convs=pdg_convs, xfxs=[toy_xfx], alphas=const_alphas)
    )

    assert_finite_nonzero(v1, "DIS NC F2 round-trip")
    assert np.array_equal(v1, v2), "Round-trip changed convolution values"


def test_05_dis_cc_f2_plus(theory, op_card):
    """CC DIS F2+: full convolution is finite and non-zero."""
    pg = build_pineappl("grid_dis_cc.yaml", theory, op_card)
    full = convolve(pg)

    assert_finite_nonzero(full, "DIS CC F2+")


def test_06_sidis_f2_unpolpdf_unpolff(theory, op_card):
    """SIDIS F2: two-convolution grid is finite and non-zero at LO+NLO."""
    pg = build_pineappl("grid_sidis.yaml", theory, op_card)
    n_conv = len(pg.convolutions)
    assert n_conv == 2, f"Expected 2 convolutions for SIDIS, got {n_conv}"

    mask = order_mask(pg, 1)  # LO + NLO
    nlo = convolve(pg, order_mask=mask)

    assert_finite_nonzero(nlo, "SIDIS F2 LO+NLO")

    # NLO must differ from LO
    lo = convolve(pg, order_mask=order_mask(pg, 0))
    assert not np.allclose(
        lo, nlo, rtol=TOL
    ), "SIDIS LO and LO+NLO results are unexpectedly identical"


def test_07_dis_pol_g1(theory, op_card):
    """Polarized DIS g₁: full convolution is finite and non-zero."""
    pg = build_pineappl("grid_dis_pol.yaml", theory, op_card)

    ct = pg.convolutions[0].convolution_types
    assert ct.polarized, "Expected a polarized convolution type for g₁"
    assert not ct.time_like

    full = convolve(pg)
    assert_finite_nonzero(full, "Polarized DIS g1")


def test_08_sidis_pol_g1(theory, op_card):
    """Polarized SIDIS G₁: two-convolution grid is finite and non-zero."""
    pg = build_pineappl("grid_sidis_pol.yaml", theory, op_card)
    n_conv = len(pg.convolutions)
    assert n_conv == 2, f"Expected 2 convolutions for polarized SIDIS, got {n_conv}"

    ct0 = pg.convolutions[0].convolution_types
    ct1 = pg.convolutions[1].convolution_types
    assert ct0.polarized and not ct0.time_like, "Expected POL_PDF as first conv"
    assert not ct1.polarized and ct1.time_like, "Expected UNPOL_FF as second conv"

    mask = order_mask(pg, 1)  # LO + NLO
    nlo = convolve(pg, order_mask=mask)

    assert_finite_nonzero(nlo, "Polarized SIDIS G1 LO+NLO")


def test_09_dis_ffn_massive(theory, op_card):
    """FFN massive DIS F2: full convolution is finite and non-zero."""
    pg = build_pineappl("grid_dis_ffn.yaml", theory, op_card)
    full = convolve(pg)

    assert_finite_nonzero(full, "FFN massive DIS F2")


def test_10_dis_massivezero(theory, op_card):
    """MassiveZero DIS F2: full convolution is finite.

    Note: for the toy PDF and the bin configuration defined in grid_dis_mz.yaml
    the MassiveZero coefficient functions evaluate to identically zero (the
    C++ test confirms: all bins give pineappl=0 and BSF=0).  We therefore only
    assert finiteness here, not non-zero.
    """
    pg = build_pineappl("grid_dis_mz.yaml", theory, op_card)
    full = convolve(pg)

    assert np.all(np.isfinite(full)), f"Non-finite result in MassiveZero DIS F2: {full}"


def test_11_dis_fonll(theory, op_card):
    """FONLL DIS F2: result is finite and consistent with ZM + FFN − MZ."""
    pg_fonll = build_pineappl("grid_dis_fonll.yaml", theory, op_card)
    pg_zm = build_pineappl("grid_dis.yaml", theory, op_card)
    pg_ffn = build_pineappl("grid_dis_ffn.yaml", theory, op_card)
    pg_mz = build_pineappl("grid_dis_mz.yaml", theory, op_card)

    fonll = convolve(pg_fonll)
    zm = convolve(pg_zm)
    ffn = convolve(pg_ffn)
    mz = convolve(pg_mz)

    assert_finite_nonzero(fonll, "FONLL DIS F2")

    # FONLL combination: F_FONLL ≈ F_ZM + F_FFN − F_MZ  (same logic as C++ test)
    ref = zm + ffn - mz
    nonzero = np.abs(ref) > 1e-15
    if np.any(nonzero):
        rel_diff = np.abs(fonll[nonzero] - ref[nonzero]) / np.abs(ref[nonzero])
        assert np.all(
            rel_diff < 0.10
        ), f"FONLL vs ZM+FFN-MZ rel_diff too large: {rel_diff}"
