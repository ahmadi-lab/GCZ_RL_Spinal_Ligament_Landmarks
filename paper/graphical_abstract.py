"""Graphical abstract (MedIA: >= 1328 px wide, readable at 13 cm).

Panel A: calibrated candidate regions on held-out vertebra S03 L4 (zones from S01 + S02).
Panel B: GCZ landmarks and their GCZ-RL corrections on the same vertebra.
Panel C: primary leave-one-subject-out results (macro MAE per held-out subject).
All values are read from the saved outputs; nothing is re-run.

Usage: python paper/graphical_abstract.py [out_stem]
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyvista as pv

ROOT = Path(__file__).resolve().parents[1]
OUT = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "paper/outputs/graphical_abstract"
sys.path.insert(0, str(ROOT / "code" / "ligament"))
import cases  # noqa: E402
import gcz  # noqa: E402
from data_io import LIGAMENTS, RAW_ROOT, read_vtk_polydata  # noqa: E402

PAT, LV, CASE = "Verse581", "L4", "S03_L4"
COL = {"ALL": "#0072B2", "PLL": "#56B4E9", "LF": "#009E73", "ISL": "#E69F00",
       "SSL": "#D55E00", "ITL": "#CC79A7", "CL": "#7B3294"}
VIEW, UP = (0.55, 1.0, 0.45), (0, 0, 1)   # LPS: from left-posterior-superior

v, f = read_vtk_polydata(RAW_ROOT / PAT / f"{PAT}_{LV}.vtk")
mesh = pv.PolyData(v, np.hstack([np.full((len(f), 1), 3), f]).ravel())
train = {(s, l): {"verts": cases.mesh(s, l, "vtk"), "gt": cases.landmarks(s, l)}
         for s, ls in {"Verse005": ["L1", "L2", "L3", "L4", "L5"], "Verse561": ["L1", "L2", "L3", "L4"]}.items()
         for l in ls}
thr = gcz.build_thresholds(train)
out = ROOT / "part2_patient_loso_14" / "outputs"
def load(fn):
    d = pd.read_csv(out / fn)
    d = d[d.case == CASE]
    return {lig: g.sort_values("point_index")[["x", "y", "z"]].to_numpy(float) for lig, g in d.groupby("ligament")}
p_gcz, p_rl = load("predictions_gcz.csv"), load("predictions_gczrl.csv")


def plotter():
    pl = pv.Plotter(off_screen=True, window_size=(1000, 1000))
    pl.set_background("white")
    pl.add_mesh(mesh, color="#ece7dd", smooth_shading=True, specular=0.1)
    return pl


def finish(pl):
    pl.view_vector(VIEW, viewup=UP)
    pl.camera.zoom(1.35)
    img = pl.screenshot(return_img=True)
    pl.close()
    return img


pl = plotter()
for lig in ["ALL", "PLL", "LF", "ISL", "CL", "SSL", "ITL"]:
    z = gcz.candidate_vertices(v, lig, thr, third_fallback=False)
    if len(z) > 4000:
        z = z[np.random.default_rng(0).choice(len(z), 4000, replace=False)]
    pl.add_points(z, color=COL[lig], point_size=4 if lig not in ("SSL", "ITL") else 7,
                  render_points_as_spheres=True, opacity=0.55 if lig not in ("SSL", "ITL") else 1.0)
img_a = finish(pl)

pl = plotter()
for lig in LIGAMENTS:
    a, b = p_gcz.get(lig), p_rl.get(lig)
    if a is None or b is None:
        continue
    pl.add_points(a, color="#9a9a9a", point_size=13, render_points_as_spheres=True)
    pl.add_points(b, color=COL[lig], point_size=17, render_points_as_spheres=True)
    for i in range(min(len(a), len(b))):
        if np.linalg.norm(a[i] - b[i]) > 0.3:
            pl.add_lines(np.vstack([a[i], b[i]]), color="#111111", width=5)
img_b = finish(pl)

subj = pd.read_csv(out / "table_P2_4_per_subject.csv")
fig = plt.figure(figsize=(13.0, 4.6), dpi=200)
gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, 1.05], wspace=0.08)
for k, (img, title) in enumerate([(img_a, "A  Calibrated admissible regions"),
                                  (img_b, "B  Bounded residual correction")]):
    ax = fig.add_subplot(gs[k])
    ax.imshow(img[60:-40, 60:-60])
    ax.set_axis_off()
    ax.set_title(title, loc="left", fontsize=13, fontweight="bold")
ax = fig.add_subplot(gs[2])
x = np.arange(3)
cols = {"GCZ": "#9a9a9a", "SLD": "#56B4E9", "GCZ-RL": "#D55E00"}
for i, m in enumerate(["GCZ", "SLD", "GCZ-RL"]):
    vals = [float(subj[(subj.method == m) & (subj.held_out_subject == s)].MAE_mm.iloc[0]) for s in ["S01", "S02", "S03"]]
    ax.bar(x + (i - 1) * 0.27, vals, width=0.25, color=cols[m], label=m)
ax.set_xticks(x, ["S01", "S02", "S03"], fontsize=12)
ax.set_ylabel("Macro mean error (mm)", fontsize=12)
ax.tick_params(axis="y", labelsize=11)
ax.set_ylim(0, 5.6)
ax.set_title("C  Leave-one-subject-out", loc="left", fontsize=13, fontweight="bold")
ax.spines[["top", "right"]].set_visible(False)
ax.yaxis.grid(True, color="#e4e3df", linewidth=0.6)
ax.set_axisbelow(True)
ax.legend(frameon=False, fontsize=11, loc="upper right", ncol=3, bbox_to_anchor=(1.0, 0.86))
ax.text(0.02, 0.985, "GCZ 2.88 → GCZ-RL 2.15 mm\nSLD 2.22 mm", transform=ax.transAxes, va="top",
        fontsize=11.5)
handles = [plt.Line2D([], [], marker="o", ls="", color=COL[l], markersize=9, label=l) for l in LIGAMENTS]
handles.append(plt.Line2D([], [], marker="o", ls="", color="#9a9a9a", markersize=8, label="GCZ point"))
fig.legend(handles=handles, loc="lower left", ncol=8, frameon=False, fontsize=11, bbox_to_anchor=(0.01, -0.01),
           handletextpad=0.2, columnspacing=0.8)
fig.subplots_adjust(left=0.01, right=0.99, top=0.9, bottom=0.12)
OUT.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(OUT.with_suffix(".png"), dpi=200)
fig.savefig(OUT.with_suffix(".pdf"))
print("saved", OUT)
