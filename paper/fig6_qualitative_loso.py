"""Qualitative figure for the primary (patient-level LOSO) evaluation.

One held-out vertebra, three orthogonal views (superior/axial, left lateral/sagittal,
posterior/coronal), expert landmarks and the predictions of GCZ, GCZ-RL and SLD.
Every point shown is read from the saved prediction files; nothing is re-run.

Usage: python paper/fig6_qualitative_loso.py [out_png] [patient] [level]
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
OUT = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "paper/outputs/fig6_qualitative_loso_S03_L4.png"
PAT = sys.argv[2] if len(sys.argv) > 2 else "Verse581"
LV = sys.argv[3] if len(sys.argv) > 3 else "L4"
sys.path.insert(0, str(ROOT / "code" / "ligament"))
import cases  # noqa: E402
from data_io import LIGAMENTS, RAW_ROOT, SUBJECT_CODE, read_vtk_polydata  # noqa: E402

COL = {"ALL": "#0072B2", "PLL": "#56B4E9", "LF": "#009E73", "ISL": "#E69F00",
       "SSL": "#D55E00", "ITL": "#CC79A7", "CL": "#7B3294"}
METHODS = [("GCZ", "predictions_gcz.csv"), ("GCZ-RL", "predictions_gczrl.csv"), ("SLD", "predictions_sld.csv")]
# LPS frame: +x left, +y posterior, +z superior
VIEWS = [("Superior (axial)", (0, 0, 1), (0, -1, 0)),
         ("Left lateral (sagittal)", (1, 0, 0), (0, 0, 1)),
         ("Posterior (coronal)", (0, 1, 0), (0, 0, 1))]

v, f = read_vtk_polydata(RAW_ROOT / PAT / f"{PAT}_{LV}.vtk")
mesh = pv.PolyData(v, np.hstack([np.full((len(f), 1), 3), f]).ravel())
gt = cases.landmarks(PAT, LV)
case = f"{SUBJECT_CODE[PAT]}_{LV}"
preds = {}
for name, fn in METHODS:
    df = pd.read_csv(ROOT / "part2_patient_loso_14" / "outputs" / fn)
    df = df[df.case == case]
    preds[name] = {lig: g[["x", "y", "z"]].to_numpy(float) for lig, g in df.groupby("ligament")}


def render(points_by_lig, expert, view, up):
    pl = pv.Plotter(off_screen=True, window_size=(900, 900))
    pl.set_background("white")
    pl.add_mesh(mesh, color="#e9e4da", smooth_shading=True, specular=0.15)
    for lig in LIGAMENTS:
        g = expert.get(lig)
        if g is not None and len(g):
            pl.add_points(np.atleast_2d(g), color="#4d4d4d", point_size=6, render_points_as_spheres=True)
    for lig in LIGAMENTS:
        p = points_by_lig.get(lig)
        if p is not None and len(p):
            pl.add_points(np.atleast_2d(p), color=COL[lig], point_size=17, render_points_as_spheres=True)
    pl.view_vector(view, viewup=up)
    pl.camera.zoom(1.25)
    img = pl.screenshot(return_img=True)
    pl.close()
    return img


fig, axes = plt.subplots(len(METHODS), len(VIEWS), figsize=(9.0, 9.2))
for r, (name, _) in enumerate(METHODS):
    for c, (vname, vec, up) in enumerate(VIEWS):
        ax = axes[r, c]
        ax.imshow(render(preds[name], gt, vec, up))
        ax.set_xticks([])
        ax.set_yticks([])
        for s in ax.spines.values():
            s.set_visible(False)
        if r == 0:
            ax.set_title(vname, fontsize=10)
        if c == 0:
            ax.set_ylabel(name, fontsize=11, fontweight="bold")
handles = [plt.Line2D([], [], marker="o", ls="", color=COL[l], markersize=8, label=l) for l in LIGAMENTS]
handles.append(plt.Line2D([], [], marker="o", ls="", color="#4d4d4d", markersize=4, label="Expert"))
fig.legend(handles=handles, loc="lower center", ncol=8, frameon=False, fontsize=9, handletextpad=0.2,
           columnspacing=0.9)
fig.tight_layout(rect=(0, 0.04, 1, 1))
OUT.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(OUT, dpi=300)
print("saved", OUT, case)
