"""Plot combined PineAPFEL vs APFEL++ accuracy for DIS, SIA, and SIDIS.

Run from project root:
    python3 plotting/plot_accuracy_combined.py
"""

import pathlib
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict

MP_STYLE = pathlib.Path("benchmarks/mplstyle.mplstyle")
if MP_STYLE.exists():
    plt.style.use(str(MP_STYLE))

ORDER_COLORS = {"LO": "C0", "NLO": "C1", "NNLO": "C2"}
ORDER_LS = {"LO": "-", "NLO": "--", "NNLO": "-."}


def parse_dat(path):
    """Parse accuracy_*.dat → dict keyed by (Q2_str, order)."""
    data = defaultdict(lambda: {"x": [], "rd": []})
    if not pathlib.Path(path).exists():
        return None
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            cols = line.split()
            q2, x, order = float(cols[0]), float(cols[1]), cols[2]
            rd = float(cols[5])
            key = (q2, order)
            data[key]["x"].append(x)
            data[key]["rd"].append(rd)
    for key in data:
        data[key]["x"] = np.asarray(data[key]["x"])
        data[key]["rd"] = np.asarray(data[key]["rd"])
    return data


def make_combined_figure(output_file):
    processes = [
        (
            "plotting/accuracy_dis.dat",
            r"$\mathrm{DIS}$ $F_2^{\mathrm{NC}}$ $\mathrm{ZM}$",
            "$x$",
        ),
        (
            "plotting/accuracy_sia.dat",
            r"$\mathrm{SIA}$ $F_2^{\mathrm{NC}}$ $\mathrm{ZM}$",
            "$z$",
        ),
        (
            "plotting/accuracy_sidis.dat",
            r"$\mathrm{SIDIS}$ $F_T^{\mathrm{NC}}$ $\mathrm{ZM}$",
            "$x$",
        ),
    ]

    q2_targets = [10.0, 100.0]

    nrows = len(processes)
    ncols = len(q2_targets)

    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=(5.0 * ncols, 3.5 * nrows),
        squeeze=False,
        sharex="row",
        sharey="all",
    )

    for row_idx, (dat_file, title_prefix, xlabel) in enumerate(processes):
        data = parse_dat(dat_file)
        if data is None:
            print(f"  {dat_file} not found — skipping row {row_idx}")
            continue

        for col_idx, Q2 in enumerate(q2_targets):
            ax = axes[row_idx][col_idx]
            any_plotted = False

            # Find closest Q2 in data
            available_q2s = sorted(list({k[0] for k in data.keys()}))
            if not available_q2s:
                continue

            # Find the Q2 in data that is closest to our target
            Q2_actual = min(available_q2s, key=lambda x: abs(x - Q2))
            if abs(Q2_actual - Q2) > 0.1:  # Allow some tolerance
                print(f"Q2={Q2} not found in {dat_file}. Closest {Q2_actual}")

            for order in ("LO", "NLO", "NNLO"):
                key = (Q2_actual, order)
                if key not in data:
                    continue
                d = data[key]
                xv = d["x"]
                rd = d["rd"]
                srt = np.argsort(xv)
                ax.plot(
                    xv[srt],
                    rd[srt],
                    color=ORDER_COLORS[order],
                    ls=ORDER_LS[order],
                    lw=1.5,
                    label=order,
                )
                any_plotted = True

            if any_plotted:
                # Tolerance reference lines
                for tol, ls, lbl in [
                    (1e-3, "k--", r"$10^{-3}$"),
                    (1e-4, "k:", r"$10^{-4}$"),
                ]:
                    ax.axhline(
                        tol,
                        ls=ls.replace("k", ""),
                        color="gray",
                        lw=0.8,
                        alpha=0.6,
                    )

            ax.set_xscale("log")
            ax.set_yscale("log")

            if row_idx == nrows - 1:
                ax.set_xlabel(xlabel, fontsize=12)

            if col_idx == 0:
                ax.set_ylabel(
                    r"$|\Delta F / F|$",
                    fontsize=12,
                )

            Q2_label = f"$Q^2 = {Q2_actual:.0f}~\\mathrm{{GeV}}^2$"
            ax.set_title(f"{title_prefix} ({Q2_label})", fontsize=11)

            if row_idx == 0 and col_idx == ncols - 1:
                ax.legend(fontsize=9, ncol=3, loc="upper right")

            # ax.grid(True, which="both", ls=":", lw=0.4, alpha=0.5)
            ax.set_xlim(1e-3, 1.0)
            ax.set_ylim(1e-15, 0)

    fig.tight_layout()
    fig.savefig(output_file, bbox_inches="tight", dpi=450)
    print(f"Saved {output_file}")
    plt.close(fig)


if __name__ == "__main__":
    print("=== PineAPFEL combined accuracy plot ===")
    make_combined_figure("plotting/accuracy_combined.pdf")
