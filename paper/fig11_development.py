"""Fig. 11 of the manuscript: development analysis (values of Table 9), drawn from
part1_dissertation_57/outputs/table_4_4_heldout_results.csv in the design of the original figure
(per-ligament MAE with the 5 mm line, and Pass@5, for the original GCZ, GCZ-RL, and SLD).

Run part1_dissertation_57/run_part1.py first.
Usage: python paper/fig11_development.py
Output: paper/outputs/figures/fig11_development_results.{pdf,svg,png}
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import figstyle as fs  # noqa: E402

REPO = Path(__file__).resolve().parents[1]
TABLE = REPO / "part1_dissertation_57" / "outputs" / "table_4_4_heldout_results.csv"
OUT = REPO / "paper" / "outputs" / "figures"
METHODS = [("GCZ", "Original GCZ", "#3d76c4", "#0b2a8a"),
           ("GCZ-RL", "GCZ-RL (Extra Trees)", "#e3241b", "#e3241b"),
           ("SLD", "SLD (same-dataset reproduction)", "#b4d4f2", "#2a6fbf")]
CATS = ["ALL", "PLL", "LF", "ISL", "SSL", "ITL", "CL", "Mean"]


def main():
    fs.setup()
    plt = fs.plt
    t = pd.read_csv(TABLE)
    val = {(r.method, r.ligament): (r.MAE_mm, r.Pass5_percent) for r in t.itertuples()}
    assert all((m, c) in val for m, *_ in METHODS for c in CATS)
    x = np.arange(len(CATS))
    fig, (a, b) = plt.subplots(2, 1, figsize=(7.0, 4.9), sharex=True, gridspec_kw={"hspace": 0.12})
    fig.subplots_adjust(left=0.075, right=0.99, top=0.985, bottom=0.11)
    for ax, k, fmt, ylim, step in ((a, 0, "{:.3f}", 10, 1), (b, 1, "{:.1f}%", 125, 25)):
        for i, (m, _, face, text) in enumerate(METHODS):
            xs = x + (i - 1) * 0.285
            h = [val[(m, c)][k] for c in CATS]
            ax.bar(xs, h, width=0.25, color=face, edgecolor="none", zorder=2)
            for xi, hi in zip(xs, h):
                ax.text(xi, hi + ylim * 0.012, fmt.format(hi), ha="center", va="bottom", fontsize=5.2, color=text, zorder=3,
                        bbox=dict(fc="white", ec="none", pad=0.15, alpha=0.85))
        ax.set_ylim(0, ylim)
        ax.set_yticks(np.arange(0, ylim + 0.1, step))
        ax.yaxis.grid(True, color="#e6e6e6", linewidth=0.6, linestyle=(0, (4, 3)))
        ax.set_axisbelow(True)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
    a.axhline(5, color=fs.INK, lw=0.9, ls=(0, (5, 3)), zorder=1)
    a.text(len(CATS) - 0.55, 5.15, "5 mm", ha="right", va="bottom", fontsize=7.5, style="italic")
    a.set_ylabel("MAE (mm)")
    b.set_ylabel("Pass rate @ 5 mm (%)")
    b.set_xticks(x, CATS, fontweight="bold")
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    a.legend(handles=[Patch(fc=face, label=lab) for _, lab, face, _ in METHODS]
             + [Line2D([], [], color=fs.INK, lw=0.9, ls=(0, (5, 3)), label="5 mm acceptance threshold")],
             loc="upper center", frameon=True, framealpha=1.0, edgecolor="#cccccc", fontsize=7)
    fig.text(0.53, 0.012, "Pass@5 mm means the percentage of Hungarian-matched landmark distances that are ≤ 5 mm.",
             ha="center", va="bottom", fontsize=6.8, style="italic")
    files = fs.save(fig, OUT, "fig11_development_results")
    print("Fig. 11:", ", ".join(p.name for p in files))


if __name__ == "__main__":
    main()
