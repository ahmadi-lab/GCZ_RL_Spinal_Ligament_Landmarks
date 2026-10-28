"""Figures of the revised manuscript, drawn from the saved data of the code release.

Every anatomical panel uses the real meshes and the saved landmark coordinates; every plotted
value is read from the outputs of the analysis scripts (run part2_patient_loso_14/run_part2.py,
paper/make_paper_numbers.py, paper/revision_analyses.py, and paper/icp_seed_sensitivity.py
first). The representative vertebra is S03 L4 of the fold that holds out S03, the vertebra of
Fig. 6; its regions are calibrated on S01 and S02, and its atlas is S01 L1, as in the primary run.

  fig1_pipeline            pipeline on real data (calibration, regions, ICP transfer, GCZ, GCZ-RL)
  fig3_calibration         calibration bands and candidate regions
  fig4_residual_correction GCZ initialization, GCZ-RL refinement, and expert reference
  fig5_patient_level       per held-out subject and per vertebra
  fig6_qualitative         GCZ, GCZ-RL, SLD on S03 L4 in three views
  fig7_components          component analyses
  fig8_icp_seeds           ten-seed ICP sensitivity
  fig9_development         development analysis (Table 9)
  figS3_residual_vectors   predicted and true residual vectors (Supplementary S10)

Writes paper/outputs/figures/ (PDF, SVG for quantitative figures, PNG) and figure_audit.csv.
Usage: python paper/figures_main.py [figure ...]
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "code" / "ligament"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import cases  # noqa: E402
import gcz  # noqa: E402
import figstyle as fs  # noqa: E402
from data_io import RAW_ROOT, SUBJECT_CODE, read_vtk_polydata  # noqa: E402
from metrics import table_from_csv  # noqa: E402

plt = fs.plt
OUT = REPO / "paper" / "outputs" / "figures"
P2 = REPO / "part2_patient_loso_14" / "outputs"
REV = REPO / "paper" / "outputs" / "revision"
SEEDS = REPO / "paper" / "outputs" / "icp_seed_sensitivity"
TARGET, ATLAS = ("Verse581", "L4"), ("Verse005", "L1")
TRAIN = [("Verse005", lv) for lv in ("L1", "L2", "L3", "L4", "L5")] + [("Verse561", lv) for lv in ("L1", "L2", "L3", "L4")]
TARGET_SEED = 200 + 1000 + 3        # full-ICP seed of S03 L4 in run_part2.gcz_predict (seed0 1000, index 3)
CORRECTED = ["ISL", "SSL", "ITL", "CL"]   # classes with blend weight > 0 in the fold holding out S03
AUDIT = []


def audit(fig, files, sources, note=""):
    AUDIT.append({"figure": fig, "files": "; ".join(p.name for p in files),
                  "sources": "; ".join(str(s.relative_to(REPO)) if isinstance(s, Path) else s for s in sources),
                  "note": note})


def mesh(key):
    v, f = read_vtk_polydata(RAW_ROOT / key[0] / f"{key[0]}_{key[1]}.vtk")
    return np.asarray(v, float), f


def code(key):
    return f"{SUBJECT_CODE[key[0]]} {key[1]}"


def preds(name):
    return table_from_csv(P2 / f"predictions_{name}.csv")


def calibration():
    train = {k: {"verts": cases.mesh(*k, "vtk"), "gt": cases.landmarks(*k)} for k in TRAIN}
    return train, gcz.build_thresholds(train)


def region_masks(v, thr, third_fallback=False):
    kd = cKDTree(v)
    out = {}
    for lig in fs.LIGAMENTS:
        zone = gcz.candidate_vertices(v, lig, thr, third_fallback=third_fallback)
        d, idx = kd.query(np.asarray(zone, float))
        m = np.zeros(len(v), bool)
        m[idx[d < 1e-6]] = True
        out[lig] = m
    return out


def class_spheres(points_by_lig, size=24, ligs=fs.LIGAMENTS):
    return [{"kind": "spheres", "points": points_by_lig[l], "color": fs.CLASS[l], "size": size}
            for l in ligs if points_by_lig.get(l) is not None and len(points_by_lig[l])]


def class_cubes(points_by_lig, edge=2.6, ligs=fs.LIGAMENTS, opacity=1.0):
    return [{"kind": "cubes", "points": points_by_lig[l], "color": fs.CLASS[l], "edge": edge, "opacity": opacity}
            for l in ligs if points_by_lig.get(l) is not None and len(points_by_lig[l])]


def expert_layer(gt, ligs=fs.LIGAMENTS, size=9):
    pts = [np.atleast_2d(gt[l]) for l in ligs if gt.get(l) is not None and len(gt[l])]
    return {"kind": "spheres", "points": np.vstack(pts), "color": fs.EXPERT, "size": size}


def correction_pairs(gcz_pts, rl_pts, ligs):
    pairs = []
    for l in ligs:
        a, b = gcz_pts.get(l), rl_pts.get(l)
        if a is None or b is None:
            continue
        a, b = np.atleast_2d(a), np.atleast_2d(b)
        pairs += list(zip(a, b))           # the two files store each class in the same point order
    return pairs


def class_legend(fig, extra=(), y=0.0, ncol=10):
    items = [(l, fs.CLASS[l], "o", 6) for l in fs.LIGAMENTS] + list(extra)
    fig.legend(handles=fs.legend_handles(items), loc="lower center", ncol=ncol, bbox_to_anchor=(0.5, y),
               handletextpad=0.2, columnspacing=0.9)


# ------------------------------------------------------------------------------ Fig. 1
def fig1():
    v, f = mesh(TARGET)
    av, af = mesh(ATLAS)
    train, thr = calibration()
    regions = region_masks(v, thr)
    atlas_gt = cases.landmarks(*ATLAS)
    r, t = gcz.icp(av, v, n_sub=5000, max_iter=150, tol=1e-6, seed=TARGET_SEED)
    av_al = gcz.apply(r, t, av)
    transferred = {l: gcz.apply(r, t, atlas_gt[l]) for l in fs.LIGAMENTS if atlas_gt.get(l) is not None}
    ref = next(k for k in TRAIN if k[1] == "L3")
    fb = gcz.build_fallback_template(train, ref, dict(n_sub=5000, max_iter=100, tol=1e-6, seed=42))
    recomputed = gcz.predict_from_alignment(v, r, t, av, atlas_gt, thr, fb)
    g, rl = preds("gcz")[TARGET], preds("gczrl")[TARGET]
    dev = max(max(cKDTree(np.atleast_2d(recomputed[l])).query(np.atleast_2d(g[l]))[0].max(),
                  cKDTree(np.atleast_2d(g[l])).query(np.atleast_2d(recomputed[l]))[0].max()) for l in g)
    gt = cases.landmarks(*TARGET)
    focal = v.mean(0)
    view, zoom = "superior", 1.45
    imgs = [
        fs.render(av, af, view, class_spheres(atlas_gt, size=16), zoom=zoom),
        fs.render(v, f, view, zoom=zoom, vertex_rgb=fs.region_colors(v, regions)),
        fs.render(v, f, view, [{"kind": "mesh", "verts": av_al, "faces": af, "color": "#5f7f9f", "opacity": 0.35}]
                  + class_spheres(transferred, size=16), zoom=zoom, focal=focal),
        fs.render(v, f, view, class_cubes(g), zoom=zoom, vertex_rgb=fs.region_colors(v, regions, mix=0.3)),
        fs.render(v, f, view, [expert_layer(gt, size=12), {"kind": "lines", "pairs": correction_pairs(g, rl, fs.LIGAMENTS),
                                                  "color": "#3a3a3a", "radius": 0.3}] + class_spheres(rl, size=40),
                  zoom=zoom)]
    titles = ["Calibration annotations\n(training vertebra S01 L1)", "Calibrated candidate regions\n(target S03 L4)",
              "Atlas transfer by rigid ICP\n(S01 L1 aligned to S03 L4)", "GCZ landmarks snapped\ninto the regions",
              "GCZ-RL refinement and\nexpert reference"]
    fig, axes = plt.subplots(1, 5, figsize=(7.2, 2.05))
    for ax, img, letter, title in zip(axes, imgs, "abcde", titles):
        fs.show(ax, img)
        ax.set_title(f"({letter}) {title}", fontsize=6.8, pad=3, loc="left")
    class_legend(fig, [("Expert reference", fs.EXPERT, "o", 3), ("GCZ landmark", fs.MUTED, "s", 5),
                       ("GCZ-RL landmark", fs.MUTED, "o", 6)], y=-0.05, ncol=10)
    fig.subplots_adjust(left=0.005, right=0.995, top=0.82, bottom=0.12, wspace=0.04)
    files = fs.save(fig, OUT, "fig1_pipeline", vector=False)
    audit("Fig. 1", files, [P2 / "predictions_gcz.csv", P2 / "predictions_gczrl.csv",
                            "dataset/raw (VerSe_Ligaments_GT meshes and annotations of S01, S02, S03 L4)"],
          f"ICP recomputed with the run seed; snapped landmarks equal the saved GCZ points (max deviation {dev:.2e} mm)")
    return dev


# ------------------------------------------------------------------------------ Fig. 3
def fig3():
    v, f = mesh(TARGET)
    train, thr = calibration()
    regions = region_masks(v, thr)
    vals = {l: [] for l in fs.LIGAMENTS}
    for d in train.values():
        for l in fs.LIGAMENTS:
            g = (d["gt"] or {}).get(l)
            if g is not None and len(g):
                vals[l].append(gcz.threshold_features(gcz.aabb_norm(g, d["verts"]))["y_n"])
    vals = {l: np.concatenate(x) for l, x in vals.items() if x}
    fig = plt.figure(figsize=(7.2, 4.3))
    gs = fig.add_gridspec(2, 3, height_ratios=[1, 1.15], hspace=0.38, wspace=0.12)
    ax = fig.add_subplot(gs[0, 0:2])
    b = thr["PLL"]["y_n"]
    ax.axvspan(b["lo"], b["hi"], color=fs.GRID, zorder=0)
    ax.hist(vals["PLL"], bins=np.arange(0, 1.0001, 0.02), color=fs.CLASS["PLL"], edgecolor="white", linewidth=0.4)
    ax.axvline(b["mu"], color=fs.INK, linewidth=0.8)
    for x in (b["lo"], b["hi"]):
        ax.axvline(x, color=fs.INK, linewidth=0.6, linestyle=(0, (3, 2)))
    ax.text(b["mu"], ax.get_ylim()[1] * 0.97, r"$\mu_L$", ha="left", va="top", fontsize=7)
    ax.text(b["hi"], ax.get_ylim()[1] * 0.97, r" $\mu_L+2\sigma_L$", ha="left", va="top", fontsize=7)
    ax.text(b["lo"], ax.get_ylim()[1] * 0.97, r"$\mu_L-2\sigma_L$ ", ha="right", va="top", fontsize=7)
    ax.set_xlim(0.15, 0.70)
    ax.set_xlabel(r"Normalized antero-posterior coordinate $y_n$ of the expert PLL points")
    ax.set_ylabel("Points")
    fs.panel(ax, "a", f"Calibration of the PLL band ({len(vals['PLL'])} points, 9 training vertebrae)")
    ax2 = fig.add_subplot(gs[0, 2])
    for i, l in enumerate(fs.LIGAMENTS[::-1]):
        bb = thr[l]["y_n"]
        ax2.plot([bb["lo"], bb["hi"]], [i, i], color=fs.CLASS[l], linewidth=2.2, solid_capstyle="butt")
        ax2.plot([bb["mu"]], [i], marker="|", color=fs.INK, markersize=7)
    ax2.set_yticks(range(len(fs.LIGAMENTS)), fs.LIGAMENTS[::-1])
    ax2.set_xlim(0, 1)
    ax2.set_xlabel(r"$y_n$ band ($\mu_L \pm 2\sigma_L$)")
    fs.panel(ax2, "b", "Bands of all classes")
    for j, view in enumerate(("superior", "left", "posterior")):
        a = fig.add_subplot(gs[1, j])
        fs.show(a, fs.render(v, f, view, zoom=1.35, vertex_rgb=fs.region_colors(v, regions)))
        a.set_title(f"({'cde'[j]}) {fs.VIEW_TITLE[view]}", fontsize=8, fontweight="bold", loc="left", pad=2)
    class_legend(fig, [("Outside all regions", fs.BONE, "s", 6)], y=-0.01, ncol=8)
    fig.subplots_adjust(left=0.07, right=0.99, top=0.95, bottom=0.08)
    files = fs.save(fig, OUT, "fig3_calibration", vector=False)
    audit("Fig. 3", files, ["dataset/raw (annotations of S01 L1-L5, S02 L1-L4; mesh of S03 L4)",
                            "code/ligament/gcz.py build_thresholds, candidate_vertices"],
          "Regions are the projection regions of the fold holding out S03 (Eq. 4, k_L = 2)")


# ------------------------------------------------------------------------------ Fig. 4
def fig4():
    v, f = mesh(TARGET)
    _, thr = calibration()
    regions = region_masks(v, thr)
    corr_regions = {l: regions[l] for l in CORRECTED}
    g, rl = preds("gcz")[TARGET], preds("gczrl")[TARGET]
    gt = cases.landmarks(*TARGET)
    focal = np.vstack([np.atleast_2d(g[l]) for l in CORRECTED]).mean(0)
    layers = [expert_layer(gt, CORRECTED, size=11),
              {"kind": "lines", "pairs": correction_pairs(g, rl, CORRECTED), "color": "#202020", "radius": 0.35}]
    layers += class_spheres(rl, size=36, ligs=CORRECTED) + class_cubes(g, edge=2.2, ligs=CORRECTED, opacity=0.45)
    rgb = fs.region_colors(v, corr_regions, mix=0.22)
    m = pd.read_csv(P2 / "per_vertebra_per_ligament_metrics.csv")
    m = m[m.case == f"{TARGET[0]}_{TARGET[1]}"].pivot_table(index="ligament", columns="method", values="mae_mm")
    fig = plt.figure(figsize=(7.2, 2.75))
    gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, 0.82], wspace=0.24)
    for j, view in enumerate(("posterior", "left")):
        a = fig.add_subplot(gs[0, j])
        fs.show(a, fs.render(v, f, view, layers, zoom=1.9, focal=focal, vertex_rgb=rgb))
        a.set_title(f"({'ab'[j]}) {fs.VIEW_TITLE[view]}", fontsize=8, fontweight="bold", loc="left", pad=2)
    ax = fig.add_subplot(gs[0, 2])
    y = np.arange(len(fs.LIGAMENTS))[::-1]
    for yi, l in zip(y, fs.LIGAMENTS):
        a_, b_ = m.loc[l, "GCZ"], m.loc[l, "GCZ-RL"]
        ax.plot([a_, b_], [yi, yi], color=fs.NEUTRAL, linewidth=1.0, zorder=1)
        ax.scatter(a_, yi, marker=fs.MARKER["GCZ"], s=18, color=fs.METHOD["GCZ"], zorder=2, linewidths=0)
        ax.scatter(b_, yi, marker=fs.MARKER["GCZ-RL"], s=20, color=fs.METHOD["GCZ-RL"], zorder=3, linewidths=0)
    ax.set_yticks(y, fs.LIGAMENTS)
    ax.set_xlabel("MAE on S03 L4 (mm)")
    ax.xaxis.grid(True, color=fs.GRID, linewidth=0.5)
    ax.set_axisbelow(True)
    fs.panel(ax, "c", "Error by class")
    ax.legend(handles=fs.legend_handles([("GCZ", fs.METHOD["GCZ"], "s", 4.5), ("GCZ-RL", fs.METHOD["GCZ-RL"], "o", 5)]),
              loc="lower right")
    items = [(l, fs.CLASS[l], "o", 6) for l in CORRECTED] + [
        ("Expert reference", fs.EXPERT, "o", 3), ("GCZ landmark", fs.MUTED, "s", 5), ("GCZ-RL landmark", fs.MUTED, "o", 6)]
    fig.legend(handles=fs.legend_handles(items), loc="lower center", ncol=7, bbox_to_anchor=(0.37, -0.02),
               handletextpad=0.2, columnspacing=0.9)
    fig.subplots_adjust(left=0.005, right=0.99, top=0.9, bottom=0.12)
    files = fs.save(fig, OUT, "fig4_residual_correction", vector=False)
    audit("Fig. 4", files, [P2 / "predictions_gcz.csv", P2 / "predictions_gczrl.csv",
                            P2 / "per_vertebra_per_ligament_metrics.csv", "dataset/raw (S03 L4 mesh and annotations)"],
          "Classes whose blend weight was above zero in the fold holding out S03")


# ------------------------------------------------------------------------------ Fig. 5
def fig5():
    subj = pd.read_csv(P2 / "table_P2_4_per_subject.csv")
    vert = pd.read_csv(REPO / "paper" / "outputs" / "primary_per_vertebra_macro_mae.csv")
    vert["label"] = [f"{SUBJECT_CODE[c.split('_')[0]]}\n{c.split('_')[1]}" for c in vert.case]
    fig, (a, b) = plt.subplots(1, 2, figsize=(7.2, 2.35), gridspec_kw={"width_ratios": [1, 3.2], "wspace": 0.22})
    order = ["GCZ", "SLD", "GCZ-RL"]
    for i, s in enumerate(["S01", "S02", "S03"]):
        vals = {mth: float(subj[(subj.method == mth) & (subj.held_out_subject == s)].MAE_mm.iloc[0]) for mth in order}
        a.plot([i, i], [vals["GCZ"], vals["GCZ-RL"]], color=fs.NEUTRAL, linewidth=1.0, zorder=1)
        for k, mth in enumerate(order):
            a.scatter(i + (k - 1) * 0.12, vals[mth], marker=fs.MARKER[mth], s=24, color=fs.METHOD[mth], zorder=3,
                      linewidths=0, label=mth if i == 0 else None)
    a.set_xticks(range(3), ["S01", "S02", "S03"])
    a.set_xlim(-0.5, 2.5)
    a.set_ylabel("Macro-averaged MAE (mm)")
    a.set_xlabel("Held-out subject")
    fs.ygrid(a)
    fs.panel(a, "a", "By held-out subject")
    x = np.arange(len(vert))
    for xi, a_, b_ in zip(x, vert["GCZ"], vert["GCZ-RL"]):
        b.plot([xi, xi], [a_, b_], color=fs.NEUTRAL, linewidth=1.0, zorder=1)
    b.scatter(x - 0.18, vert.GCZ, marker="s", s=18, color=fs.METHOD["GCZ"], linewidths=0, zorder=3)
    b.scatter(x + 0.18, vert.SLD, marker="^", s=18, color=fs.METHOD["SLD"], linewidths=0, zorder=2)
    b.scatter(x, vert["GCZ-RL"], marker="o", s=20, color=fs.METHOD["GCZ-RL"], linewidths=0, zorder=3)
    b.set_xticks(x, vert.label, fontsize=6.3)
    for s in (4.5, 8.5):
        b.axvline(s, color=fs.GRID, linewidth=0.8)
    b.set_xlim(-0.6, len(vert) - 0.4)
    fs.ygrid(b)
    fs.panel(b, "b", "By held-out vertebra (GCZ-RL lower than GCZ for 13 of 14)")
    a.legend(loc="upper right", handletextpad=0.1, borderaxespad=0.1)
    files = fs.save(fig, OUT, "fig5_patient_level")
    audit("Fig. 5", files, [P2 / "table_P2_4_per_subject.csv", REPO / "paper/outputs/primary_per_vertebra_macro_mae.csv"])


# ------------------------------------------------------------------------------ Fig. 6
def fig6():
    v, f = mesh(TARGET)
    gt = cases.landmarks(*TARGET)
    methods = [("GCZ", preds("gcz")[TARGET]), ("GCZ-RL", preds("gczrl")[TARGET]), ("SLD", preds("sld")[TARGET])]
    fig, axes = plt.subplots(3, 3, figsize=(6.4, 6.3))
    for r_, (name, p) in enumerate(methods):
        for c, view in enumerate(("superior", "left", "posterior")):
            ax = axes[r_, c]
            fs.show(ax, fs.render(v, f, view, [expert_layer(gt, size=9)] + class_spheres(p, size=26), zoom=1.3))
            if r_ == 0:
                ax.set_title(fs.VIEW_TITLE[view], fontsize=8)
            if c == 0:
                ax.set_ylabel(name, fontsize=8.5, fontweight="bold")
    class_legend(fig, [("Expert reference", fs.EXPERT, "o", 3)], y=0.0, ncol=8)
    fig.subplots_adjust(left=0.04, right=0.995, top=0.96, bottom=0.05, wspace=0.02, hspace=0.04)
    files = fs.save(fig, OUT, "fig6_qualitative", vector=False)
    audit("Fig. 6", files, [P2 / "predictions_gcz.csv", P2 / "predictions_gczrl.csv", P2 / "predictions_sld.csv",
                            "dataset/raw (S03 L4)"], "Same points as paper/check_fig6.py")


# ------------------------------------------------------------------------------ Fig. 7
ABLATION_LABELS = ["GCZ (geometry only)", "GCZ with projection, no residual", "GCZ with raw residual",
                   "GCZ with bounded residual", r"GCZ with bounded residual and projection ($\alpha_L=1$)",
                   r"Complete GCZ-RL (nested $\alpha_L$)"]
# Drawing styles of Fig. 7: A, monochrome; B, the method colours of the other figures (used in the paper).
FIG7_STYLE = {"A": {"GCZ": ("#8a8a8a", "s"), "GCZ-RL": (fs.INK, "o"), "other": ("#4d4d4d", "o")},
              "B": {"GCZ": (fs.METHOD["GCZ"], "s"), "GCZ-RL": (fs.METHOD["GCZ-RL"], "o"), "other": ("#4d4d4d", "o")}}
CLASS_LABEL_OFFSET = {"CL": (-5, 0, "right")}


def axis_break(ax, which="y"):
    """Mark an axis that does not start at zero with two short slanted strokes near its origin."""
    kw = dict(transform=ax.transAxes, color=fs.INK, lw=0.7, clip_on=False, zorder=5)
    gap = dict(transform=ax.transAxes, color="white", lw=1.6, clip_on=False, zorder=4)
    if which == "y":
        ax.plot([0, 0], [0.023, 0.047], **gap)
        for d in (-0.012, 0.012):
            ax.plot([-0.022, 0.022], [0.035 + d - 0.012, 0.035 + d + 0.012], **kw)
    else:
        ax.plot([0.023, 0.047], [0, 0], **gap)
        for d in (-0.012, 0.012):
            ax.plot([0.035 + d - 0.012, 0.035 + d + 0.012], [-0.03, 0.03], **kw)


def fig7(variant="B", stem="fig7_components"):
    abl = pd.read_csv(REV / "table_ablation_components.csv")
    comp = pd.read_csv(REV / "table_region_compactness.csv")
    perm = pd.read_csv(REV / "feature_importance_permutation.csv")
    atl = pd.read_csv(REV / "table_atlas_icp.csv")
    st = FIG7_STYLE[variant]
    fig = plt.figure(figsize=(7.2, 5.3))
    gs = fig.add_gridspec(2, 2, hspace=0.66, wspace=0.62)
    fig.subplots_adjust(left=0.335, right=0.985, top=0.93, bottom=0.085)

    # (a) Cleveland dot plot of the ablation, in methodological order
    ax = fig.add_subplot(gs[0, 0])
    vals = abl.macro_mae_mm.to_numpy()
    n = len(vals)
    keys = ["GCZ"] + ["other"] * (n - 2) + ["GCZ-RL"]
    ax.yaxis.grid(True, color=fs.GRID, linewidth=0.6, linestyle=(0, (1, 1.5)))
    ax.set_axisbelow(True)
    for i, (v, k) in enumerate(zip(vals, keys)):
        col, mk = st[k]
        ax.plot(v, i, marker=mk, ms=5.4 if k != "other" else 4.6, color=col, ls="", zorder=3)
        ax.text(1.04, i, f"{v:.3f}", transform=ax.get_yaxis_transform(), va="center", ha="left", fontsize=7,
                color=fs.INK, fontweight="bold" if k != "other" else "normal")
    ax.text(1.04, -0.95, "MAE (mm)", transform=ax.get_yaxis_transform(), va="bottom", ha="left", fontsize=6.3,
            color=fs.MUTED)
    ax.set_yticks(range(n), ABLATION_LABELS)
    for t, k in zip(ax.get_yticklabels(), keys):
        if k != "other":
            t.set_fontweight("bold")
    ax.set_ylim(n - 0.5, -0.6)
    ax.set_xlim(1.9, 3.0)
    ax.set_xticks(np.arange(2.0, 3.01, 0.2))
    ax.set_xlabel("Macro-averaged MAE (mm)")
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    axis_break(ax, "x")

    # (b) containment against compactness, one uniform marker per class with direct labels
    ax = fig.add_subplot(gs[0, 1])
    for r in comp.itertuples():
        x, y = r.projection_area_fraction_median, r.inside_selection_pct
        ax.plot(x, y, "o", ms=4.4, color=st["other"][0], ls="", zorder=3)
        dx, dy, ha = CLASS_LABEL_OFFSET.get(r.ligament, (5, 0, "left"))
        ax.annotate(r.ligament, (x, y), xytext=(dx, dy), textcoords="offset points", ha=ha, va="center", fontsize=7)
    ax.set_xlim(-1, 30)
    ax.set_ylim(80, 100)
    ax.set_yticks(range(80, 101, 5))
    ax.set_xlabel("Share of the vertebral surface in the region (%)")
    ax.set_ylabel("Reference landmarks inside (%)")
    fs.ygrid(ax)
    axis_break(ax, "y")
    fs.panel(ax, "b", "Region coverage and compactness", x=-0.2)

    # (c) grouped permutation importance (unchanged content)
    ax = fig.add_subplot(gs[1, 0])
    g = perm.groupby("group").increase_mm
    order = g.mean().sort_values().index
    yy = np.arange(len(order))
    ax.barh(yy, g.mean()[order], height=0.5, color=fs.NEUTRAL)
    ax.errorbar(g.mean()[order], yy, xerr=[g.mean()[order] - g.min()[order], g.max()[order] - g.mean()[order]],
                fmt="none", ecolor=fs.INK, elinewidth=0.6, capsize=1.5)
    ax.set_yticks(yy, [o.replace(" (level, class)", "") for o in order])
    ax.set_xlabel("Increase in held-out residual error (mm)")
    ax.xaxis.grid(True, color=fs.GRID, linewidth=0.5)
    ax.set_axisbelow(True)
    p = ax.get_position()  # shorter labels than (a): widen the axes to the left
    ax.set_position([0.2, p.y0, p.x1 - 0.2, p.height])

    # (d) atlas residual against per-vertebra error; legend above the axes, outside the data
    ax = fig.add_subplot(gs[1, 1])
    for r in atl.itertuples():
        ax.plot([r.icp_residual_mm] * 2, [r.GCZ_mm, r.GCZRL_mm], color=fs.GRID, linewidth=0.9, zorder=1)
    for same in (True, False):
        d = atl[atl.same_level == same]
        for col, mth in (("GCZ_mm", "GCZ"), ("GCZRL_mm", "GCZ-RL")):
            c, mk = st[mth]
            ax.scatter(d.icp_residual_mm, d[col], marker=mk, s=20, facecolors=c if same else "white", edgecolors=c,
                       linewidths=0.8, zorder=3)
    ax.set_xlabel("Quick-ICP residual of the selected atlas (mm)")
    ax.set_ylabel("Per-vertebra MAE (mm)")
    fs.ygrid(ax)
    from matplotlib.lines import Line2D
    ax.legend(handles=[Line2D([], [], ls="", marker=st["GCZ"][1], color=st["GCZ"][0], label="GCZ", markersize=4.5),
                       Line2D([], [], ls="", marker=st["GCZ-RL"][1], color=st["GCZ-RL"][0], label="GCZ-RL",
                              markersize=4.5),
                       Line2D([], [], ls="", marker="o", markerfacecolor="white", markeredgecolor=fs.MUTED,
                              label="Cross-level atlas (open markers)", markersize=4.5)],
              loc="lower left", bbox_to_anchor=(-0.02, 1.0), ncol=3, frameon=False, handletextpad=0.15,
              columnspacing=0.9, borderaxespad=0.15, fontsize=6.6)
    ax.set_title("(d) Atlas residual and landmark error", loc="left", fontsize=8, fontweight="bold", x=-0.2, pad=17)

    # titles of the left panels start at the figure margin, above their long category labels
    for ax_, title in ((fig.axes[0], "(a) Component ablation"), (fig.axes[2], "(c) Grouped permutation importance")):
        p = ax_.get_position()
        ax_.set_title(title, loc="left", fontsize=8, fontweight="bold", x=(0.012 - p.x0) / p.width, pad=6)
    files = fs.save(fig, OUT, stem)
    audit("Fig. 7" if stem == "fig7_components" else f"Fig. 7 alternative {variant}", files,
          [REV / "table_ablation_components.csv", REV / "table_region_compactness.csv",
           REV / "feature_importance_permutation.csv", REV / "table_atlas_icp.csv"],
          "" if stem == "fig7_components" else "Design alternative for comparison; not used in the manuscript")


def fig7_alternatives():
    """Alternative A (monochrome) and B (two method colours) for comparison; B is used as Fig. 7."""
    for v in ("A", "B"):
        fig7(v, f"fig7_components_alt{v}")


# ------------------------------------------------------------------------------ Fig. 8
def fig8():
    s = pd.read_csv(SEEDS / "summary.csv")
    stats = json.loads((SEEDS / "summary.json").read_text())
    rep = json.loads((REPO / "paper" / "outputs" / "paper_numbers.json").read_text())
    fig, ax = plt.subplots(figsize=(3.5, 2.5))
    xs = [-1.4] + list(s.offset)
    gz = [rep["GCZ"]["macro_mae"]] + list(s.GCZ_mm)
    rl = [rep["GCZ-RL"]["macro_mae"]] + list(s.GCZRL_mm)
    for x, a, b in zip(xs, gz, rl):
        ax.plot([x, x], [a, b], color=fs.NEUTRAL, linewidth=1.0, zorder=1)
    ax.scatter(xs[1:], gz[1:], marker="s", s=16, color=fs.METHOD["GCZ"], linewidths=0, zorder=3, label="GCZ")
    ax.scatter(xs[1:], rl[1:], marker="o", s=18, color=fs.METHOD["GCZ-RL"], linewidths=0, zorder=3, label="GCZ-RL")
    ax.scatter(xs[:1], gz[:1], marker="s", s=16, facecolors="white", edgecolors=fs.METHOD["GCZ"], linewidths=0.9, zorder=3)
    ax.scatter(xs[:1], rl[:1], marker="o", s=18, facecolors="white", edgecolors=fs.METHOD["GCZ-RL"], linewidths=0.9,
               zorder=3)
    sld = rep["SLD"]["macro_mae"]
    ax.axhline(sld, color=fs.METHOD["SLD"], linewidth=0.8, linestyle=(0, (4, 2)))
    ax.text(0.5, sld + 0.03, "SLD (no ICP)", color=fs.METHOD["SLD"], fontsize=6.3, ha="center", va="bottom")
    ax.set_xticks(xs, ["Stored"] + [str(int(k)) for k in s.offset], fontsize=6.5)
    ax.set_xlabel("ICP seed offset (alignments recomputed)")
    ax.set_ylabel("Macro-averaged MAE (mm)")
    ax.set_ylim(1.9, 4.1)
    fs.ygrid(ax)
    red = stats["reduction_pct"]
    ax.text(0.02, 0.98, f"GCZ-RL < GCZ in {stats['runs']}/{stats['runs']} runs and all subjects\n"
                        f"Mean reduction {red['mean']:.1f}% (range {red['min']:.1f}–{red['max']:.1f}%)",
            transform=ax.transAxes, va="top", fontsize=6.5)
    ax.legend(loc="lower right", bbox_to_anchor=(1.0, 1.0), ncol=2, handletextpad=0.1, fontsize=6.5)
    files = fs.save(fig, OUT, "fig8_icp_seeds")
    audit("Fig. 8", files, [SEEDS / "summary.csv", SEEDS / "summary.json", REPO / "paper/outputs/paper_numbers.json"],
          "Computational sensitivity analysis on the same three subjects; not independent validation")


# ------------------------------------------------------------------------------ Fig. 9
def fig9():
    t = pd.read_csv(REPO / "part1_dissertation_57" / "outputs" / "table_4_4_heldout_results.csv")
    cats = fs.LIGAMENTS + ["Mean"]
    fig, axes = plt.subplots(2, 1, figsize=(3.5, 3.5), sharex=True, gridspec_kw={"hspace": 0.25})
    for ax, col, lab in ((axes[0], "MAE_mm", "MAE (mm)"), (axes[1], "Pass5_percent", "Pass@5 (%)")):
        for k, mth in enumerate(["GCZ", "SLD", "GCZ-RL"]):
            d = t[t.method == mth].set_index("ligament").reindex(cats)
            ax.scatter(np.arange(len(cats)) + (k - 1) * 0.22, d[col], marker=fs.MARKER[mth], s=16,
                       color=fs.METHOD[mth], linewidths=0, label={"GCZ": "Original GCZ", "SLD": "SLD (context only)",
                                                                   "GCZ-RL": "GCZ-RL"}[mth])
        ax.set_ylabel(lab)
        fs.ygrid(ax)
        ax.axvline(len(cats) - 1.5, color=fs.GRID, linewidth=0.8)
    axes[1].set_xticks(range(len(cats)), cats)
    axes[0].legend(loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=3, handletextpad=0.1, columnspacing=0.6,
                   fontsize=6.3)
    files = fs.save(fig, OUT, "fig9_development")
    audit("Fig. 9", files, [REPO / "part1_dissertation_57/outputs/table_4_4_heldout_results.csv"])


# ------------------------------------------------------------------------------ Fig. S3
def figS3():
    v, f = mesh(TARGET)
    pts = pd.read_csv(REV / "residual_vectors_heldout.csv")
    pts = pts[(pts.case == "S03_L4") & (pts.alpha > 0)]
    gt = cases.landmarks(*TARGET)
    p0 = pts[["gcz_x", "gcz_y", "gcz_z"]].to_numpy()
    true_end = p0 + pts[["true_x", "true_y", "true_z"]].to_numpy()
    pred_end = p0 + pts[["pred_x", "pred_y", "pred_z"]].to_numpy()
    layers = [expert_layer(gt, CORRECTED, size=11),
              {"kind": "lines", "pairs": list(zip(p0, true_end)), "color": "#2b2b2b", "radius": 0.25}]
    for l, d in pts.groupby("ligament"):
        q0 = d[["gcz_x", "gcz_y", "gcz_z"]].to_numpy()
        layers.append({"kind": "lines", "pairs": list(zip(q0, q0 + d[["pred_x", "pred_y", "pred_z"]].to_numpy())),
                       "color": fs.CLASS[l], "radius": 0.2})
        layers.append({"kind": "cubes", "points": q0, "color": fs.CLASS[l], "edge": 2.4})
    focal = p0.mean(0)
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.0))
    for ax, view, letter in zip(axes, ("posterior", "left"), "ab"):
        fs.show(ax, fs.render(v, f, view, layers, zoom=1.9, focal=focal))
        ax.set_title(f"({letter}) {fs.VIEW_TITLE[view]}", fontsize=8, fontweight="bold", loc="left", pad=2)
    from matplotlib.lines import Line2D
    items = [Line2D([], [], color=fs.CLASS[l], lw=1.5, marker="s", markersize=4, label=f"{l}: GCZ point, predicted residual")
             for l in CORRECTED]
    items += [Line2D([], [], color="#2b2b2b", lw=1.5, label="True residual"),
              Line2D([], [], ls="", marker="o", color=fs.EXPERT, markersize=3, label="Expert reference")]
    fig.legend(handles=items, loc="lower center", ncol=3, bbox_to_anchor=(0.5, -0.04))
    fig.subplots_adjust(left=0.005, right=0.995, top=0.92, bottom=0.17, wspace=0.03)
    files = fs.save(fig, OUT, "figS3_residual_vectors", vector=False)
    audit("Fig. S3", files, [REV / "residual_vectors_heldout.csv", "dataset/raw (S03 L4)"],
          "Vectors at true length, before the bound and projection")


def main():
    fs.setup()
    names = sys.argv[1:] or ["fig1", "fig3", "fig4", "fig5", "fig6", "fig7", "fig7_alternatives", "fig8", "fig9",
                             "figS3"]
    for n in names:
        res = globals()[n]()
        print(n, "done" if res is None else f"done ({res})", flush=True)
    path = OUT / "figure_audit.csv"
    old = pd.read_csv(path) if path.exists() else pd.DataFrame(columns=["figure", "files", "sources", "note"])
    new = pd.DataFrame(AUDIT)
    pd.concat([old[~old.figure.isin(new.figure)], new]).sort_values("figure").to_csv(path, index=False)


if __name__ == "__main__":
    main()
