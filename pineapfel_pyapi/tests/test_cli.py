"""
CLI regression tests for pineapfel-build and pineapfel-evolve.

Each test invokes the CLI binary, convolves the resulting grid with the toy
distribution, and compares against a stored .npy reference.

Run with --generate-references once to create / refresh the .npy files:
    pytest pineapfel_pyapi/tests/test_cli.py --generate-references -v
"""

import os
import subprocess

import numpy as np

from conftest import convolve_path, runcard


REFERENCES = os.path.join(os.path.dirname(__file__), "references")


def ref_path(name: str) -> str:
    return os.path.join(REFERENCES, name)


def _run(cmd: list[str], env: dict | None = None) -> None:
    """Run a subprocess, raising on non-zero exit."""
    result = subprocess.run(cmd, capture_output=True, text=True, env=env)
    if result.returncode != 0:
        raise RuntimeError(f"Command {cmd} failed:\n{result.stdout}\n{result.stderr}")


def _check_or_generate(
    values: np.ndarray, filename: str, generate: bool, rtol: float = 1e-6
) -> None:
    path = ref_path(filename)
    if generate:
        os.makedirs(REFERENCES, exist_ok=True)
        np.save(path, values)
        return
    ref = np.load(path)
    np.testing.assert_allclose(
        values,
        ref,
        rtol=rtol,
        err_msg=f"CLI regression failed for {filename}",
    )


def test_cli_build_dis_nc(
    pineapfel_build_bin, theory_path, op_path, generate, tmp_path
):
    """pineapfel-build: DIS NC F2 grid matches reference convolution."""
    out = str(tmp_path / "dis_nc.pineappl.lz4")
    _run(
        [
            pineapfel_build_bin,
            runcard("grid_dis.yaml"),
            theory_path,
            op_path,
            "-o",
            out,
        ]
    )
    values = convolve_path(out)
    _check_or_generate(values, "cli_dis_nc_f2.npy", generate)


def test_cli_build_dis_cc(
    pineapfel_build_bin, theory_path, op_path, generate, tmp_path
):
    """pineapfel-build: DIS CC F2+ grid matches reference convolution."""
    out = str(tmp_path / "dis_cc.pineappl.lz4")
    _run(
        [
            pineapfel_build_bin,
            runcard("grid_dis_cc.yaml"),
            theory_path,
            op_path,
            "-o",
            out,
        ]
    )
    values = convolve_path(out)
    _check_or_generate(values, "cli_dis_cc_f2.npy", generate)


def test_cli_build_dis_pol(
    pineapfel_build_bin, theory_path, op_path, generate, tmp_path
):
    """pineapfel-build: polarized DIS g1 grid matches reference convolution."""
    out = str(tmp_path / "dis_pol.pineappl.lz4")
    _run(
        [
            pineapfel_build_bin,
            runcard("grid_dis_pol.yaml"),
            theory_path,
            op_path,
            "-o",
            out,
        ]
    )
    values = convolve_path(out)
    _check_or_generate(values, "cli_dis_pol_g1.npy", generate)


# ---------------------------------------------------------------------------
# pineapfel-evolve tests
# ---------------------------------------------------------------------------


def test_cli_evolve_dis_nc(
    pineapfel_build_bin, pineapfel_evolve_bin, theory_path, op_path, generate, tmp_path
):
    """pineapfel-build + pineapfel-evolve: DIS NC FK table matches reference."""
    grid_out = str(tmp_path / "dis_nc.pineappl.lz4")
    fk_out = str(tmp_path / "dis_nc.fk.pineappl.lz4")

    _run(
        [
            pineapfel_build_bin,
            runcard("grid_dis.yaml"),
            theory_path,
            op_path,
            "-o",
            grid_out,
        ]
    )
    _run([pineapfel_evolve_bin, grid_out, theory_path, op_path, "-o", fk_out])

    values = convolve_path(fk_out)
    _check_or_generate(values, "cli_fk_dis_nc_f2.npy", generate)
