"""Figures of the revision analyses (reads paper/outputs/revision, written by revision_analyses.py).

  fig_component_analysis.{pdf,png}   (a) component ablation, (b) region coverage vs compactness,
                                     (c) grouped permutation importance, (d) ICP residual vs GCZ error
  fig_residual_vectors_S03_L4.png    true and predicted residual vectors on a held-out vertebra

Usage: python paper/revision_figures.py
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
REV = ROOT / "paper" / "outputs" / "revision"
sys.path.insert(0, str(ROOT / "code" / "ligament"))
import cases  # noqa: E402
from data_io import LIGAMENTS, RAW_ROOT, SUBJECT_CODE, read_vtk_polydata  # noqa: E402

COL = {"ALL": "#0072B2", "PLL": "#56B4E9", "LF": "#009E73", "ISL": "#E69F00",
       "SSL": "#D55E00", "ITL": "#CC79A7", "CL": "#7B3294"}
METHOD = {"GCZ": "#2a78d6", "GCZ-RL": "#eb6834", "SLD": "#1baf7a"}
INK, MUTED, GRID, STAGE_GRAY = "#262624", "#6b6a66", "#e4e3df", "#a9a7a1"
plt.rcParams.update({"font.size": 8, "axes.edgecolor": MUTED, "axes.labelcolor": INK,
                     "xtick.color": MUTED, "ytick.color": MUTED, "font.family": "DejaVu Sans"})


def style(ax):
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_axisbelow(True)


def components():
    abl = pd.read_csv(REV / "table_ablation_components.csv")
    comp = pd.read_csv(REV / "table_region_compactness.csv")
    perm = pd.read_csv(REV / "feature_importance_permutation.csv")
    atl = pd.read_csv(REV / "table_atlas_icp.csv")
    fig, axes = plt.subplots(2, 2, figsize=(7.4, 6.4))

    ax = axes[0, 0]
    labels = ["GCZ", "+ projection only", "+ raw residual", "+ bounded residual", "+ bounded + projection",
              "Complete GCZ-RL"]
    colors = [METHOD["GCZ"]] + [STAGE_GRAY] * 4 + [METHOD["GCZ-RL"]]
    y = np.arange(len(abl))[::-1]
    ax.barh(y, abl.macro_mae_mm, height=0.62, color=colors, edgecolor="white", linewidth=0.8)
    for yi, (_, row) in zip(y, abl.iterrows()):
        ax.text(row.macro_mae_mm + 0.05, yi, f"{row.macro_mae_mm:.3f}  ({row.admissible_vertex_pct:.0f}%)",
                va="center", ha="left", color=INK, fontsize=6.6)
    ax.set_yticks(y, labels, fontsize=7)
    ax.set_xlabel("Macro-averaged MAE (mm)")
    ax.set_xlim(0, 4.3)
    ax.xaxis.grid(True, color=GRID, linewidth=0.6)
    ax.set_title("(a) Component ablation", fontsize=8, loc="left", color=INK)
    style(ax)

    ax = axes[0, 1]
    for _, row in comp.iterrows():
        ax.scatter(row.projection_area_fraction_median, row.inside_selection_pct, s=42, color=COL[row.ligament],
                   edgecolor="white", linewidth=1.2, zorder=3)
        ax.annotate(row.ligament, (row.projection_area_fraction_median, row.inside_selection_pct),
                    xytext=(5, 3), textcoords="offset points", fontsize=7, color=INK)
    ax.set_xlabel("Region share of surface area (%, median)")
    ax.set_ylabel("Reference landmarks inside region (%)")
    ax.set_xlim(0, max(comp.projection_area_fraction_median) * 1.25 + 1)
    ax.set_ylim(80, 100)
    ax.grid(True, color=GRID, linewidth=0.6)
    ax.set_title("(b) Region coverage versus compactness", fontsize=8, loc="left", color=INK)
    style(ax)

    ax = axes[1, 0]
    g = perm.groupby("group").increase_mm
    order = g.mean().sort_values().index
    y = np.arange(len(order))
    mean, lo, hi = g.mean()[order], g.min()[order], g.max()[order]
    ax.barh(y, mean, height=0.55, color=METHOD["GCZ-RL"], edgecolor="white")
    ax.errorbar(mean, y, xerr=[mean - lo, hi - mean], fmt="none", ecolor=INK, elinewidth=0.8, capsize=2)
    ax.set_yticks(y, order, fontsize=7)
    ax.set_xlabel("Increase in held-out residual error (mm)")
    ax.xaxis.grid(True, color=GRID, linewidth=0.6)
    ax.set_title("(c) Grouped permutation importance", fontsize=8, loc="left", color=INK)
    style(ax)

    ax = axes[1, 1]
    for same, mk in ((True, "o"), (False, "s")):
        d = atl[atl.same_level == same]
        ax.scatter(d.icp_residual_mm, d.GCZ_mm, marker=mk, s=34, color=METHOD["GCZ"], edgecolor="white",
                   linewidth=1.0, zorder=3, label="GCZ, same-level atlas" if same else "GCZ, cross-level atlas")
        ax.scatter(d.icp_residual_mm, d.GCZRL_mm, marker=mk, s=34, color=METHOD["GCZ-RL"], edgecolor="white",
                   linewidth=1.0, zorder=3, label="GCZ-RL, same-level atlas" if same else "GCZ-RL, cross-level atlas")
    for _, row in atl.iterrows():
        ax.plot([row.icp_residual_mm] * 2, [row.GCZ_mm, row.GCZRL_mm], color=GRID, linewidth=1.0, zorder=1)
        if row.atlas_changed:
            ax.annotate(row.target.replace("_", " "), (row.icp_residual_mm, row.GCZ_mm), xytext=(6, 2),
                        textcoords="offset points", fontsize=6.5, color=INK)
    ax.set_xlabel("Quick-ICP residual of the selected atlas (mm)")
    ax.set_ylabel("Per-vertebra macro MAE (mm)")
    ax.grid(True, color=GRID, linewidth=0.6)
    ax.legend(frameon=False, fontsize=6.2, loc="upper center", bbox_to_anchor=(0.45, -0.19), ncol=2, handletextpad=0.2,
              columnspacing=0.8)
    ax.set_ylim(1.4, 4.7)
    ax.set_title("(d) Atlas residual versus landmark error", fontsize=8, loc="left", color=INK)
    style(ax)

    fig.tight_layout(w_pad=3.0, h_pad=1.5)
    for ext in ("pdf", "png"):
        fig.savefig(REV / f"fig_component_analysis.{ext}", dpi=300)
    plt.close(fig)


def vectors(pat="Verse581", lv="L4"):
    pts = pd.read_csv(REV / "residual_vectors_heldout.csv")
    case = f"{SUBJECT_CODE[pat]}_{lv}"
    pts = pts[(pts.case == case) & (pts.alpha > 0)]
    v, f = read_vtk_polydata(RAW_ROOT / pat / f"{pat}_{lv}.vtk")
    mesh = pv.PolyData(v, np.hstack([np.full((len(f), 1), 3), f]).ravel())
    gt = cases.landmarks(pat, lv)
    views = [("Posterior (coronal)", (0, 1, 0), (0, 0, 1)), ("Left lateral (sagittal)", (1, 0, 0), (0, 0, 1))]
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.6))
    for ax, (name, vec, up) in zip(axes, views):
        pl = pv.Plotter(off_screen=True, window_size=(1100, 900))
        pl.set_background("white")
        pl.add_mesh(mesh, color="#e9e4da", smooth_shading=True, specular=0.15, opacity=0.85)
        for lig in pts.ligament.unique():
            g = gt.get(lig)
            if g is not None and len(g):
                pl.add_points(np.atleast_2d(g), color="#4d4d4d", point_size=5, render_points_as_spheres=True)
        for lig, d in pts.groupby("ligament"):
            p0 = d[["gcz_x", "gcz_y", "gcz_z"]].to_numpy()
            t = d[["true_x", "true_y", "true_z"]].to_numpy()
            pr = d[["pred_x", "pred_y", "pred_z"]].to_numpy()
            pl.add_points(p0, color=COL[lig], point_size=14, render_points_as_spheres=True)
            pl.add_arrows(p0, t, mag=1.0, color="#262624")
            pl.add_arrows(p0, pr, mag=1.0, color=COL[lig])
        pl.view_vector(vec, viewup=up)
        pl.camera.focal_point = tuple(np.vstack([d[["gcz_x", "gcz_y", "gcz_z"]].to_numpy() for _, d in pts.groupby("ligament")]).mean(0))
        pl.camera.zoom(1.6)
        img = pl.screenshot(return_img=True)
        rows, cols = np.where(img.min(axis=2) < 250)
        ax.imshow(img[max(rows.min() - 10, 0):rows.max() + 10, max(cols.min() - 10, 0):cols.max() + 10])
        pl.close()
        ax.set_xticks([])
        ax.set_yticks([])
        for s in ax.spines.values():
            s.set_visible(False)
        ax.set_title(name, fontsize=9, color=INK)
    handles = [plt.Line2D([], [], marker="o", ls="", color=COL[l], markersize=7, label=f"{l} GCZ point / predicted residual")
               for l in pts.ligament.unique()]
    handles += [plt.Line2D([], [], color="#262624", lw=2, label="True residual (to paired reference)"),
                plt.Line2D([], [], marker="o", ls="", color="#4d4d4d", markersize=4, label="Expert reference")]
    fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False, fontsize=7)
    fig.tight_layout(rect=(0, 0.13, 1, 0.97))
    fig.savefig(REV / f"fig_residual_vectors_{case}.png", dpi=300)
    plt.close(fig)


if __name__ == "__main__":
    components()
    vectors()
    print("figures written to", REV)
