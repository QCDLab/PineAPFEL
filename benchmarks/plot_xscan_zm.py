"""Plot absolute relative differences from benchmarks/xscan_zm.dat.

Run from the project root:
    python benchmarks/plot_xscan_zm.py

input : benchmarks/xscan_zm.dat
output: benchmarks/xscan_zm.pdf
"""

import pathlib

import matplotlib.pyplot as plt
import matplotlib.colors as clr
import numpy as np

DAT_FILE = pathlib.Path("benchmarks/xscan_zm.dat")
PDF_FILE = pathlib.Path("benchmarks/xscan_zm.pdf")
MP_STYLE = pathlib.Path("benchmarks/mplstyle.mplstyle")

if MP_STYLE.exists():
    plt.style.use(str(MP_STYLE))

ORDER_MARKERS = {"LO": "o", "NLO": "s", "NNLO": "o"}
ORDERS_PER_PROC = {
    "DIS": ["LO", "NLO", "NNLO"],
    "SIA": ["LO", "NLO", "NNLO"],
    "SIDIS": ["LO", "NLO"],
}
MAP_ORDERS_LABEL = {
    "LO": r"$\mathrm{LO}$",
    "NLO": r"$\mathrm{NLO}$",
    "NNLO": r"$\mathrm{NNLO}$",
}

DATA: dict = {}

with open(DAT_FILE) as fh:
    for line in fh:
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        proc_ord = parts[0]
        sep = proc_ord.rfind("_")
        proc = proc_ord[:sep]
        order = proc_ord[sep + 1 :]
        xz = float(parts[1])
        pine = float(parts[2])
        apfel = float(parts[3])
        ratio = pine / apfel if apfel != 0.0 else 1.0

        key = (proc, order)
        if key not in DATA:
            DATA[key] = {"x": [], "ratio": []}
        DATA[key]["x"].append(xz)
        DATA[key]["ratio"].append(ratio)

for key in DATA:
    for field in DATA[key]:
        DATA[key][field] = np.array(DATA[key][field])

fig, axes = plt.subplots(
    1,
    3,
    figsize=(13.5, 3.75),
    sharey=True,
    gridspec_kw={"wspace": 0},
)

process_labels = {
    "DIS": r"$\mathrm{DIS}~F_2^{\rm NC}~\mathrm{ZMVFNS~at}~Q^2=10$",
    "SIA": r"$\mathrm{SIA}~F_2^{\rm NC}~\mathrm{ZMVFNS~at}~Q^2=10$",
    "SIDIS": r"$\mathrm{SIDIS}~F_2^{\rm NC}~\mathrm{ZMVFNS~at}~Q^2=10$",
}
x_labels = {
    "DIS": r"$x$",
    "SIA": r"$z$",
    "SIDIS": r"$x$",
}

for ax, proc in zip(axes, ("DIS", "SIA", "SIDIS")):
    for order in ORDERS_PER_PROC[proc]:
        key = (proc, order)
        if key not in DATA:
            continue
        d = DATA[key]
        ax.scatter(
            d["x"],
            d["ratio"],
            s=8,
            marker=ORDER_MARKERS[order],
            label=MAP_ORDERS_LABEL[order],
            zorder=3,
        )

    ax_band = 1e-10
    ax_refb = 5e-11
    ax.set_xscale("log")
    ax.grid(False)
    ax.axhspan(
        1 - ax_refb,
        1 + ax_refb,
        linewidth=1.0,
        facecolor=clr.to_rgba("red", alpha=0.05),
        hatch="\\",
        edgecolor=clr.to_rgba("red", alpha=1.00),
        zorder=0,
        label=r"$\delta = 10^{-11}$",
    )
    ax.set_ylim([1 - ax_band, 1 + ax_band])
    ax.set_xlabel(x_labels[proc])
    if proc == "DIS":
        ax.set_ylabel(r"$\mathrm{PineAPFEL~Grids} / \mathrm{APFEL++}$")
    ax.set_title(process_labels[proc])
    ax.legend(fontsize=9, markerscale=1.5, ncol=2)

fig.tight_layout()
fig.savefig(PDF_FILE, bbox_inches="tight", dpi=350)
print(f"Saved {PDF_FILE}")
