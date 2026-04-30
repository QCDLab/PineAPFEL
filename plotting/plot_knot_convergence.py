"""Plot Lagrange interpolation convergence from plotting/knot_convergence.dat.

Run from the project root:
    python3 plotting/plot_knot_convergence.py

Input : plotting/knot_convergence.dat
Output: plotting/knot_convergence.pdf
"""

import pathlib
import numpy as np
import matplotlib.pyplot as plt

DAT_FILE = pathlib.Path("plotting/knot_convergence.dat")
PDF_FILE = pathlib.Path("plotting/knot_convergence.pdf")
MP_STYLE = pathlib.Path("benchmarks/mplstyle.mplstyle")

if MP_STYLE.exists():
    plt.style.use(str(MP_STYLE))

# Columns: config  n_knots  n_nodes  x  order  rel_diff
# key: (config, order) → {"n_nodes": [], "n_knots": [], "x": [], "rd": []}
data: dict = {}

with open(DAT_FILE) as fh:
    for line in fh:
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        cfg, nk, nn, x, order, rd = line.split()
        key = (cfg, order)
        if key not in data:
            data[key] = {"n_knots": [], "n_nodes": [], "x": [], "rd": []}
        data[key]["n_knots"].append(int(nk))
        data[key]["n_nodes"].append(int(nn))
        data[key]["x"].append(float(x))
        data[key]["rd"].append(float(rd))

for key in data:
    for field in data[key]:
        data[key][field] = np.asarray(data[key][field])

ORDER_COLORS = {"LO": "C0", "NLO": "C1", "NNLO": "C2"}
ORDER_LABELS = {
    "LO": r"$\mathrm{LO}$",
    "NLO": r"$\mathrm{NLO}$",
    "NNLO": r"$\mathrm{NNLO}$",
}
CFG_LS = {"1sg": "-", "2sg": "--"}
CFG_LABELS = {"1sg": "1 subgrid", "2sg": "2 subgrids"}

ORDERS = ["LO", "NLO", "NNLO"]

fig, (ax_conv, ax_x) = plt.subplots(1, 2, figsize=(11.5, 4.5))

for order in ORDERS:
    for cfg in ("1sg", "2sg"):
        key = (cfg, order)
        if key not in data:
            continue
        d = data[key]
        nn_vals = np.unique(d["n_nodes"])
        max_rd = np.array([d["rd"][d["n_nodes"] == nn].max() for nn in nn_vals])

        label = None
        ax_conv.plot(
            nn_vals,
            max_rd,
            color=ORDER_COLORS[order],
            ls=CFG_LS[cfg],
            marker="o" if cfg == "1sg" else "s",
            markersize=4,
            label=label,
        )

# O(N^{-4}) reference slope anchored at the first 1sg NNLO point
key_ref = ("1sg", "NNLO")
if key_ref in data:
    d_ref = data[key_ref]
    nn_vals = np.unique(d_ref["n_nodes"])
    max_rd = np.array([d_ref["rd"][d_ref["n_nodes"] == nn].max() for nn in nn_vals])
    n0, r0 = nn_vals[0], max_rd[0]
    n_slope = np.linspace(0, 250, 100)
    r_slope = r0 * (n_slope / n0) ** (-4)
    ax_conv.plot(
        n_slope,
        r_slope,
        "k:",
        lw=1.2,
        label=r"$\mathcal{O}(N^{-4})$ (degree-3 Lagrange)",
        zorder=0,
    )

for order in ORDERS:
    ax_conv.plot(
        [],
        [],
        color=ORDER_COLORS[order],
        ls="-",
        lw=1.5,
        label=ORDER_LABELS[order],
    )
ax_conv.plot([], [], color="gray", ls="-", lw=1.5, label="1 subgrid")
ax_conv.plot([], [], color="gray", ls="--", lw=1.5, label="2 subgrids")

ax_conv.set_xscale("linear")
ax_conv.set_yscale("log")
ax_conv.set_xlabel(r"Joint-grid nodes $N$", fontsize=12)
ax_conv.set_ylabel(r"$\max_x \left| \Delta F / F \right|$", fontsize=12)
ax_conv.set_title(
    r"$\mathrm{DIS~NC}~F_2~\mathrm{ZM}$"
    r" at $Q^2 = 10~\mathrm{GeV}^2$, reference at $N = 2000$",
    fontsize=11,
)
ax_conv.set_xlim(0, 250)
ax_conv.set_ylim(1e-6, 1)
ax_conv.legend(fontsize=9, ncol=2)

SELECTED_NK = [20, 60, 160]
NK_COLORS = {20: "C0", 60: "C1", 160: "C2"}

for nk in SELECTED_NK:
    for cfg in ("1sg", "2sg"):
        key = (cfg, "NNLO")
        if key not in data:
            continue
        d = data[key]
        mask = d["n_knots"] == nk
        if not mask.any():
            continue
        xv = d["x"][mask]
        rd = d["rd"][mask]
        idx = np.argsort(xv)

        label = f"$N_k={nk}$" if cfg == "1sg" else None
        ax_x.plot(
            xv[idx],
            rd[idx],
            color=NK_COLORS[nk],
            ls=CFG_LS[cfg],
            lw=1.2,
            label=label,
        )

# Add dummy lines for linestyle legend
ax_x.plot([], [], color="gray", ls="-", lw=1.2, label="1 subgrid")
ax_x.plot([], [], color="gray", ls="--", lw=1.2, label="2 subgrids")

ax_x.set_xscale("log")
ax_x.set_yscale("log")
ax_x.set_xlabel(r"$x$", fontsize=12)
ax_x.set_ylabel(r"$\left| \Delta F / F \right|$", fontsize=12)
ax_x.set_title("NNLO Interpolation Error", fontsize=11)
ax_x.legend(fontsize=9, ncol=2)
ax_x.set_xlim(1e-4, 1)
ax_x.set_ylim(1e-7, 1)

fig.tight_layout()
fig.savefig(PDF_FILE, bbox_inches="tight", dpi=450)
print(f"Saved {PDF_FILE}")
