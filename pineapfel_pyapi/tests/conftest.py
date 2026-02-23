"""
Shared pytest fixtures and helpers for PineAPFEL Python tests.

These mirror the test infrastructure of tests/test_grid_vs_apfelxx.cpp,
adapted for the Python API: grids are built with the `pineapfel` bindings
and convolved with pineappl-py instead of the C++ pineappl_capi.
"""

import math
import os
import sys
import tempfile

import numpy as np
import pineappl
import pineappl.convolutions
import pytest
import pineapfel as pf

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))


PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
RUNCARDS = os.path.join(PROJECT_ROOT, "runcards")


def runcard(name: str) -> str:
    return os.path.join(RUNCARDS, name)


def toy_f(pid: int, x: float) -> float:
    """Returns f(x), NOT x*f(x).  Matches the C++ `toy_f`."""
    apid = abs(pid)
    if pid == 21:
        return math.pow(x, -0.1) * math.pow(1.0 - x, 5.0)
    if 1 <= apid <= 5:
        return math.pow(x, 0.5) * math.pow(1.0 - x, 3.0)
    return 0.0


def toy_xfx(pid: int, x: float, q2: float) -> float:
    """PineAPPL callback: returns x * f(x).  Matches the C++ `xfx_callback`."""
    return x * toy_f(pid, x)


def const_alphas(q2: float) -> float:
    """Constant alpha_s used to isolate grid structure from running effects."""
    return 0.118


def make_pdg_convs(pg: pineappl.grid.Grid) -> list:
    """Build the list of Conv objects required by Grid.convolve from the
    grid's own convolution metadata."""
    result = []
    for c in pg.convolutions:
        ct = pineappl.convolutions.ConvType(
            polarized=c.convolution_types.polarized,
            time_like=c.convolution_types.time_like,
        )
        result.append(pineappl.convolutions.Conv(ct, c.pid))
    return result


def convolve(pg: pineappl.grid.Grid, order_mask=None) -> np.ndarray:
    """Convolve `pg` with the toy distributions.

    Works for any number of convolution functions (DIS: 1, SIDIS: 2).
    """
    pdg_convs = make_pdg_convs(pg)
    n_conv = len(pdg_convs)
    xfxs = [toy_xfx] * n_conv

    kwargs = {}
    if order_mask is not None:
        kwargs["order_mask"] = np.asarray(order_mask, dtype=bool)

    return np.asarray(
        pg.convolve(
            pdg_convs=pdg_convs,
            xfxs=xfxs,
            alphas=const_alphas,
            **kwargs,
        )
    )


@pytest.fixture(scope="session")
def theory():
    return pf.load_theory_card(runcard("theory.yaml"))


@pytest.fixture(scope="session")
def op_card():
    return pf.load_operator_card(runcard("operator.yaml"))


def build_pineappl(grid_yaml: str, theory_card, op_card_obj) -> pineappl.grid.Grid:
    """Build a PineAPFEL grid and return it as a pineappl-py Grid object."""
    gdef = pf.load_grid_def(runcard(grid_yaml))
    pg = pf.build_grid(gdef, theory_card, op_card_obj)

    with tempfile.NamedTemporaryFile(suffix=".pineappl.lz4", delete=False) as f:
        path = f.name

    try:
        pg.write(path)
        return pineappl.grid.Grid.read(path)
    finally:
        os.unlink(path)
