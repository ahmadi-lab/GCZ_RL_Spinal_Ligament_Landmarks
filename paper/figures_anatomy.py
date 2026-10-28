"""Anatomical figures of the revised manuscript, drawn from the real mesh and the saved landmarks of
held-out vertebra S03 L4 (fold holding out S03; regions calibrated on S01 and S02; atlas S01 L1).

  fig_qualitative_insets   Fig. 7: expert reference, GCZ, and GCZ-RL with magnified insets of the
                           classes whose correction is active in this fold, and the per-class errors
                           of GCZ, GCZ-RL, and SLD on this vertebra
  fig_regions_margins      Fig. 8: expert attachment annotations, calibrated search regions (GCZ
                           masks), and GCZ-RL landmarks; ALL and PLL with close-ups of the superior
                           and inferior body margins

One encoding is used in both figures and never mixed:
  expert attachment annotation     black sphere
  GCZ landmark (initial)           blue cube
  GCZ-RL landmark (final)          red sphere
  applied correction               dark segment from the GCZ to the GCZ-RL landmark
  calibrated search region (mask)  gray surface tint (Fig. 8 only)
Markers are glyphs of fixed size in millimetres, so a close-up magnifies markers and anatomy alike.
Every view uses a parallel projection with an explicit camera; a close-up is a crop of the same
render as its overview, and its frame on the overview comes from the same camera.

Checks (paper/outputs/figures/anatomy_checks.json): the drawn landmarks equal the saved predictions,
the masks are those of gcz.candidate_vertices, and the errors printed in Fig. 7 equal
part2_patient_loso_14/outputs/per_vertebra_per_ligament_metrics.csv.

Usage: python paper/figures_anatomy.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import figures_main as fm  # noqa: E402
import margin_verification as mv  # noqa: E402
from figures_main import fs, plt  # noqa: E402

OUT = fm.OUT
KEY = fm.TARGET
EXPERT_C, GCZ_C, RL_C, SLD_C, SEG_C = "#1a1a1a", fs.METHOD["GCZ"], fs.METHOD["GCZ-RL"], fs.METHOD["SLD"], "#2b2b2b"
REGION_C = "#9c9c9c"
R_EXPERT, R_RL, E_GCZ = 0.55, 0.95, 1.7          # glyph sizes (mm)
CHECKS = {}


# ------------------------------------------------------------------------------ rendering
def render(v, f, view, layers=(), rgb=None, center=None, scale=None, size=2400, project=()):
    """Parallel-projection render with an explicit camera; returns the image and the pixel
    coordinates of the world points in ``project``."""
    import pyvista as pv
    vec, up = (np.asarray(a, float) for a in fs.VIEWS[view])
    center = 0.5 * (v.min(0) + v.max(0)) if center is None else np.asarray(center, float)
    scale = 0.55 * np.ptp(v, axis=0).max() if scale is None else scale
    mesh = pv.PolyData(v, np.hstack([np.full((len(f), 1), 3), f]).ravel())
    pl = pv.Plotter(off_screen=True, window_size=(size, size))
    pl.set_background("white")
    if rgb is not None:
        mesh.point_data["rgb"] = rgb
        pl.add_mesh(mesh, scalars="rgb", rgb=True, smooth_shading=True, specular=0.05, ambient=0.25)
    else:
        pl.add_mesh(mesh, color=fs.BONE, smooth_shading=True, specular=0.05, ambient=0.25)
    for lay in layers:
        pts = np.atleast_2d(np.asarray(lay.get("points", np.zeros((0, 3))), float))
        if lay["kind"] == "spheres" and len(pts):
            g = pv.PolyData(pts).glyph(geom=pv.Sphere(radius=lay["r"], theta_resolution=18, phi_resolution=18),
                                       scale=False, orient=False)
            pl.add_mesh(g, color=lay["color"], smooth_shading=True, ambient=0.35, specular=0.2)
        elif lay["kind"] == "cubes" and len(pts):
            e = lay["edge"]
            g = pv.PolyData(pts).glyph(geom=pv.Cube(x_length=e, y_length=e, z_length=e), scale=False, orient=False)
            pl.add_mesh(g, color=lay["color"], ambient=0.35, opacity=lay.get("opacity", 1.0))
        elif lay["kind"] == "lines":
            for a, b in lay["pairs"]:
                if np.linalg.norm(np.subtract(b, a)) > 0.05:
                    pl.add_mesh(pv.Line(a, b).tube(radius=lay.get("radius", 0.18)), color=lay["color"])
    pl.camera_position = [tuple(center + vec * 1000.0), tuple(center), tuple(up)]
    pl.camera.parallel_projection = True
    pl.camera.parallel_scale = scale
    pl.renderer.ResetCameraClippingRange()
    img = pl.screenshot(return_img=True)
    k = img.shape[0] / size
    ren, xy = pl.renderer, []
    for p in project:
        ren.SetWorldPoint(*np.asarray(p, float), 1.0)
        ren.WorldToDisplay()
        dx, dy, _ = ren.GetDisplayPoint()
        xy.append((dx * k, img.shape[0] - dy * k))
    pl.close()
    return img, np.asarray(xy, float).reshape(-1, 2), img.shape[0] / (2.0 * scale)   # px per mm


def content_box(img, pad=10):
    rows, cols = np.where(img.min(axis=2) < 250)
    return (max(cols.min() - pad, 0), max(rows.min() - pad, 0),
            min(cols.max() + pad, img.shape[1]), min(rows.max() + pad, img.shape[0]))


def box_around(xy, px_mm, pad_mm, aspect=None, min_mm=None):
    """Pixel box (x0, y0, x1, y1) around projected points, padded in millimetres."""
    x0, y0 = xy.min(0) - pad_mm * px_mm
    x1, y1 = xy.max(0) + pad_mm * px_mm
    if min_mm is not None:
        cx, cy, h = 0.5 * (x0 + x1), 0.5 * (y0 + y1), 0.5 * min_mm * px_mm
        x0, x1 = min(x0, cx - h), max(x1, cx + h)
        y0, y1 = min(y0, cy - h), max(y1, cy + h)
    if aspect is not None:   # width / height
        w, h = x1 - x0, y1 - y0
        cx, cy = 0.5 * (x0 + x1), 0.5 * (y0 + y1)
        if w / h < aspect:
            w = h * aspect
        else:
            h = w / aspect
        x0, x1, y0, y1 = cx - w / 2, cx + w / 2, cy - h / 2, cy + h / 2
    return np.array([x0, y0, x1, y1])


def crop(img, b):
    """Crop a pixel box; parts outside the render are filled with the white background."""
    x0, y0, x1, y1 = np.round(b).astype(int)
    p = max(0, -x0, -y0, x1 - img.shape[1], y1 - img.shape[0])
    if p:
        img = np.pad(img, ((p, p), (p, p), (0, 0)), constant_values=255)
        x0, y0, x1, y1 = x0 + p, y0 + p, x1 + p, y1 + p
    return img[y0:y1, x0:x1]


def frame(ax, b, origin, label=None, color=fs.INK):
    from matplotlib.patches import Rectangle
    x0, y0, x1, y1 = b
    ox, oy = origin
    ax.add_patch(Rectangle((x0 - ox, y0 - oy), x1 - x0, y1 - y0, fill=False, ec=color, lw=0.8))
    if label:
        ax.text(x0 - ox + 0.06 * (x1 - x0), y0 - oy + 0.06 * (y1 - y0), label, ha="left", va="top", fontsize=7,
                fontweight="bold", color=color)


def union(*boxes):
    b = np.array(boxes, float)
    return np.array([b[:, 0].min(), b[:, 1].min(), b[:, 2].max(), b[:, 3].max()])


def show(ax, img, title=None, below=None, boxed=False):
    ax.imshow(img, interpolation="lanczos")
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(boxed)
        s.set_linewidth(0.8)
    if title:
        ax.set_title(title, loc="left", fontsize=7.5, fontweight="bold", pad=3)
    if below:
        ax.text(0.0, -0.03, below, transform=ax.transAxes, fontsize=6.6, va="top", ha="left", color=fs.INK)


def layer(kind, pts, color, size):
    key = "edge" if kind == "cubes" else "r"
    return {"kind": kind, "points": np.atleast_2d(pts), "color": color, key: size}


def stack(d, ligs):
    pts = [np.atleast_2d(d[l]) for l in ligs if d.get(l) is not None and len(d[l])]
    return np.vstack(pts) if pts else np.zeros((0, 3))


def check_saved(name, drawn, ligs):
    """The drawn landmarks equal the rows of the saved prediction file, read here without the
    loader of the code release (maximum coordinate deviation, mm)."""
    df = pd.read_csv(fm.P2 / f"predictions_{name}.csv")
    df = df[(df.patient == KEY[0]) & (df.level == KEY[1])]
    dev = 0.0
    for l in ligs:
        rows = df[df.ligament == l].sort_values("point_index")[["x", "y", "z"]].to_numpy()
        dev = max(dev, float(np.abs(np.atleast_2d(drawn[l]) - rows).max()))
    CHECKS[f"{name}_drawn_vs_saved_max_deviation_mm"] = dev
    assert dev < 1e-9, (name, dev)


# ------------------------------------------------------------------------------ Fig. 7
def fig_qualitative_insets():
    v, f = fm.mesh(KEY)
    gt = fm.cases.landmarks(*KEY)
    g0, rl = fm.preds("gcz")[KEY], fm.preds("gczrl")[KEY]
    check_saved("gcz", g0, fs.LIGAMENTS)
    check_saved("gczrl", rl, fs.LIGAMENTS)
    met = pd.read_csv(fm.P2 / "per_vertebra_per_ligament_metrics.csv")
    met = met[met.case == f"{KEY[0]}_{KEY[1]}"].pivot(index="ligament", columns="method", values="mae_mm")
    xc = 0.5 * (v[:, 0].min() + v[:, 0].max())
    pairs = [p for p in fm.correction_pairs(g0, rl, fm.CORRECTED) if np.linalg.norm(np.subtract(*p)) > 0.3]
    lays = [layer("spheres", stack(gt, fs.LIGAMENTS), EXPERT_C, R_EXPERT),
            layer("cubes", stack(g0, fs.LIGAMENTS), GCZ_C, E_GCZ),
            {"kind": "lines", "pairs": pairs, "color": SEG_C, "radius": 0.2},
            layer("spheres", stack(rl, fs.LIGAMENTS), RL_C, R_RL)]

    # inset centres: the left ITL tip, and the CL landmark with the largest applied correction
    side = lambda p: np.atleast_2d(p)[np.atleast_2d(p)[:, 0] > xc]  # noqa: E731
    itl_pts = np.vstack([side(gt["ITL"]), side(g0["ITL"]), side(rl["ITL"])])
    d_cl = np.linalg.norm(np.atleast_2d(rl["CL"]) - np.atleast_2d(g0["CL"]), axis=1)
    k = int(np.argmax(d_cl))
    cl_pts = np.vstack([np.atleast_2d(g0["CL"])[k], np.atleast_2d(rl["CL"])[k]])
    CHECKS["fig7_cl_inset_correction_mm"] = round(float(d_cl[k]), 2)
    post, xy, pxmm = render(v, f, "posterior", lays, project=np.vstack([itl_pts, cl_pts]))
    b_itl = box_around(xy[:len(itl_pts)], pxmm, 4.0, aspect=1.0, min_mm=24)
    b_cl = box_around(xy[len(itl_pts):], pxmm, 4.0, aspect=1.0, min_mm=24)
    sp_pts = np.vstack([gt["ISL"], gt["SSL"], g0["ISL"], g0["SSL"], rl["ISL"], rl["SSL"]])
    lat, xy2, pxmm2 = render(v, f, "left", lays, project=sp_pts)
    b_sp = box_around(xy2, pxmm2, 3.0, aspect=1.0)

    fig = plt.figure(figsize=(7.2, 4.55))
    gs = fig.add_gridspec(2, 4, width_ratios=[1.6, 1, 1, 1], height_ratios=[1, 0.9], wspace=0.07, hspace=0.3)
    fig.subplots_adjust(left=0.01, right=0.99, top=0.95, bottom=0.085)
    cb = union(content_box(post, pad=20), b_itl + [-6, -6, 6, 6], b_cl + [-6, -6, 6, 6])
    ax = fig.add_subplot(gs[0, 0])
    show(ax, crop(post, cb), "(a) Posterior view, all classes")
    frame(ax, b_itl, cb[:2], "b")
    frame(ax, b_cl, cb[:2], "c")

    def mae(lig):
        return f"{lig}: GCZ {met.loc[lig, 'GCZ']:.2f} mm, GCZ-RL {met.loc[lig, 'GCZ-RL']:.2f} mm"

    show(fig.add_subplot(gs[0, 1]), crop(post, b_itl), "(b) ITL, left tip", mae("ITL"), boxed=True)
    cl_side = "left" if np.atleast_2d(rl["CL"])[k, 0] > xc else "right"   # LPS: +x is the patient's left
    show(fig.add_subplot(gs[0, 2]), crop(post, b_cl), f"(c) CL, {cl_side} side", mae("CL"), boxed=True)
    show(fig.add_subplot(gs[0, 3]), crop(lat, b_sp), "(d) ISL and SSL, lateral", mae("ISL") + "\n" + mae("SSL"),
         boxed=True)

    # (e) per-class error of the three methods on this vertebra
    ax = fig.add_subplot(gs[1, 1:])
    order = fs.LIGAMENTS[::-1]
    for i, l in enumerate(order):
        a0, a1, s = met.loc[l, "GCZ"], met.loc[l, "GCZ-RL"], met.loc[l, "SLD"]
        if abs(a0 - a1) > 1e-9:
            ax.annotate("", xy=(a1, i), xytext=(a0, i),
                        arrowprops=dict(arrowstyle="-|>", color=fs.NEUTRAL, lw=0.9, shrinkA=3.5, shrinkB=3.5,
                                        mutation_scale=7))
        ax.plot(a0, i, "s", ms=4.6, color=GCZ_C, zorder=3)
        ax.plot(s, i, "^", ms=4.8, color=SLD_C, zorder=3)
        ax.plot(a1, i, "o", ms=4.8, color=RL_C, zorder=4)
        if l in ("ALL", "PLL", "LF"):
            ax.text(9.55, i, r"$\alpha_L=0$", va="center", ha="right", fontsize=6.3, color=fs.MUTED)
    ax.set_yticks(range(len(order)), order)
    ax.set_ylim(-0.6, len(order) - 0.4)
    ax.set_xlim(0, 9.6)
    ax.set_xlabel("MAE on S03 L4 (mm)")
    ax.xaxis.grid(True, color=fs.GRID, linewidth=0.5)
    ax.set_axisbelow(True)
    ax.tick_params(axis="y", length=0)
    mean = met.mean()
    ax.set_title(f"(e) Error by class on this vertebra (macro MAE: GCZ {mean['GCZ']:.3f}, GCZ-RL {mean['GCZ-RL']:.3f}, "
                 f"SLD {mean['SLD']:.3f} mm)", loc="left", fontsize=7.5, fontweight="bold", pad=3)
    CHECKS["fig7_macro_mae_S03_L4"] = {k: round(float(mean[k]), 3) for k in ("GCZ", "GCZ-RL", "SLD")}

    # bottom-left: what the markers are
    ax = fig.add_subplot(gs[1, 0])
    ax.set_axis_off()
    from matplotlib.lines import Line2D
    h = [Line2D([], [], ls="", marker="o", color=EXPERT_C, ms=3.2, label="Expert attachment annotation"),
         Line2D([], [], ls="", marker="s", color=GCZ_C, ms=4.6, label="GCZ landmark (atlas transfer)"),
         Line2D([], [], ls="", marker="o", color=RL_C, ms=4.8, label="GCZ-RL landmark (corrected)"),
         Line2D([], [], color=SEG_C, lw=1.0, label="Applied correction (GCZ to GCZ-RL)"),
         Line2D([], [], ls="", marker="^", color=SLD_C, ms=4.8, label="SLD (panel e only)")]
    ax.legend(handles=h, loc="center left", frameon=False, fontsize=6.8, handletextpad=0.5, labelspacing=0.9)
    files = fs.save(fig, OUT, "fig_qualitative_insets", vector=False)
    fm.audit("Fig. 7 (new)", files, [fm.P2 / "predictions_gcz.csv", fm.P2 / "predictions_gczrl.csv",
                                     fm.P2 / "per_vertebra_per_ligament_metrics.csv", "raw mesh and annotations of S03 L4"])
    return files


# ------------------------------------------------------------------------------ Fig. 8
def fig_regions_margins():
    v, f = fm.mesh(KEY)
    gt = fm.cases.landmarks(*KEY)
    rl = fm.preds("gczrl")[KEY]
    check_saved("gczrl", rl, fs.LIGAMENTS)
    _, thr = fm.calibration()
    masks = fm.region_masks(v, thr)
    ref = mv.region_mask(v, "ALL", mv.thresholds("Verse581"))
    CHECKS["fig8_masks_equal_margin_verification"] = bool(np.array_equal(ref, masks["ALL"]))
    walls = mv.body_walls(v)
    grey = (np.array(fs.matplotlib.colors.to_rgb(REGION_C)) * 255).astype(np.uint8)
    bone = (np.array(fs.matplotlib.colors.to_rgb(fs.BONE)) * 255).astype(np.uint8)

    def tint(lig):
        rgb = np.tile(bone, (len(v), 1))
        rgb[masks[lig]] = grey
        return rgb

    def letter(ax, s):
        ax.text(0.025, 0.97, s, transform=ax.transAxes, ha="left", va="top", fontsize=7.5, fontweight="bold",
                bbox=dict(fc="white", ec="none", pad=0.6, alpha=0.85))

    exp_l = lambda l: layer("spheres", gt[l], EXPERT_C, R_EXPERT)  # noqa: E731
    rl_l = lambda l: layer("spheres", rl[l], RL_C, R_RL)  # noqa: E731
    letters = iter("abcdefghijklmnop")

    fig = plt.figure(figsize=(7.2, 7.0))
    outer = fig.add_gridspec(2, 1, height_ratios=[1.0, 1.12], hspace=0.13)
    fig.subplots_adjust(left=0.075, right=0.995, top=0.925, bottom=0.005)

    # ALL (anterior view) and PLL (posterior view of the body, posterior arch removed for display):
    # search region with expert points, and close-ups of the superior and inferior margins at one scale
    top = outer[0].subgridspec(2, 3, width_ratios=[1.0, 1.25, 1.25], wspace=0.05, hspace=0.12)
    STRIP_W, STRIP_H = 46.0, 21.0                 # close-up size, mm
    fb = mv.body_only(v, f, walls["y_post"] + 0.5)
    for r, (lig, view, faces) in enumerate((("ALL", "anterior", f), ("PLL", "posterior", fb))):
        pp = mv.parts(*KEY, lig)
        body = v[np.unique(faces)]
        center = 0.5 * (body.min(0) + body.max(0))
        scale = 0.5 * max(np.ptp(body[:, 0]), np.ptp(body[:, 2])) * 1.1
        img_r, xy, pxmm = render(v, faces, view, [exp_l(lig)], rgb=tint(lig), center=center, scale=scale,
                                 project=np.vstack([pp["S"].mean(0), pp["I"].mean(0)]))
        img_p, _, _ = render(v, faces, view, [exp_l(lig), rl_l(lig)], center=center, scale=scale)
        boxes = [np.r_[c - [STRIP_W / 2 * pxmm, STRIP_H / 2 * pxmm], c + [STRIP_W / 2 * pxmm, STRIP_H / 2 * pxmm]]
                 for c in xy]
        cbx = union(content_box(img_r, pad=15), *boxes)
        ax = fig.add_subplot(top[r, 0])
        show(ax, crop(img_r, cbx))
        letter(ax, next(letters))
        for b in boxes:
            frame(ax, b, cbx[:2])
        ax.text(-0.04, 0.5, f"{lig}\n({'anterior' if lig == 'ALL' else 'posterior'})", transform=ax.transAxes,
                rotation=90, ha="right", va="center", fontsize=7.2, fontweight="bold")
        for c, b in enumerate(boxes):
            a = fig.add_subplot(top[r, c + 1])
            show(a, crop(img_p, b), boxed=True)
            letter(a, next(letters))
        if r == 0:
            for c, t in enumerate(["Search region and expert points", "Superior margin (close-up)",
                                   "Inferior margin (close-up)"]):
                fig.axes[c].set_title(t, loc="left", fontsize=7.3, fontweight="bold", pad=3)

    # LF, ISL, SSL, ITL, CL: one camera per class; small regions are seen end-on
    bot = outer[1].subgridspec(2, 5, wspace=0.05, hspace=0.06)
    xc = 0.5 * (v[:, 0].min() + v[:, 0].max())
    left = lambda p: np.atleast_2d(p)[np.atleast_2d(p)[:, 0] > xc]  # noqa: E731  (LPS: +x is the left)
    spec = {"LF": ("superior", np.atleast_2d(gt["LF"]), "LF, superior view"),
            "ISL": ("left", np.vstack([gt["ISL"], gt["SSL"]]), "ISL, left lateral"),
            "SSL": ("posterior", np.atleast_2d(gt["SSL"]), "SSL, tip (posterior)"),
            "ITL": ("left", left(gt["ITL"]), "ITL, left tip (lateral)"),
            "CL": ("posterior", left(gt["CL"]), "CL, left side (posterior)")}
    for c, lig in enumerate(["LF", "ISL", "SSL", "ITL", "CL"]):
        view, focus, title = spec[lig]
        img_r, xy, pxmm = render(v, f, view, [exp_l(lig)], rgb=tint(lig), project=focus)
        img_p, _, _ = render(v, f, view, [exp_l(lig), rl_l(lig)])
        box = box_around(xy, pxmm, 3.0, aspect=1.0, min_mm=24 if lig in ("SSL", "ITL") else None)
        for r, img in enumerate((img_r, img_p)):
            ax = fig.add_subplot(bot[r, c])
            show(ax, crop(img, box), boxed=True)
            letter(ax, ("ghijk" if r == 0 else "lmnop")[c])
            if r == 0:
                ax.set_title(title, loc="left", fontsize=7.0, fontweight="bold", pad=3)
            if c == 0:
                ax.text(-0.04, 0.5, ["Search region and\nexpert points", "GCZ-RL landmarks\nand expert points"][r],
                        transform=ax.transAxes, rotation=90, ha="right", va="center", fontsize=7.0)
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch, Rectangle
    h = [Line2D([], [], ls="", marker="o", color=EXPERT_C, ms=3.4, label="Expert attachment annotation"),
         Line2D([], [], ls="", marker="o", color=RL_C, ms=4.8, label="GCZ-RL landmark (final)"),
         Patch(fc=REGION_C, ec="none", label="Calibrated search region (GCZ mask)"),
         Rectangle((0, 0), 1, 1, fill=False, ec=fs.INK, lw=0.8, label="Close-up frame")]
    fig.legend(handles=h, loc="upper center", ncol=4, frameon=False, bbox_to_anchor=(0.53, 0.999), fontsize=6.9)
    CHECKS["fig8_rl_inside_region"] = {l: bool(np.isin(np.round(np.atleast_2d(rl[l]), 4),
                                                       np.round(v[masks[l]], 4)).all(axis=1).all())
                                       for l in fs.LIGAMENTS}
    files = fs.save(fig, OUT, "fig_regions_margins", vector=False)
    fm.audit("Fig. 8 (new)", files, [fm.P2 / "predictions_gczrl.csv", "raw mesh and annotations of S03 L4",
                                     "gcz.build_thresholds / candidate_vertices (fold holding out S03)"])
    return files


def main():
    fs.setup()
    names = sys.argv[1:] or ["fig_qualitative_insets", "fig_regions_margins"]
    for n in names:
        print(n, [p.name for p in globals()[n]()], flush=True)
    path = OUT / "anatomy_checks.json"
    old = json.loads(path.read_text()) if path.exists() else {}
    old.update(CHECKS)
    path.write_text(json.dumps(old, indent=2))
    print(json.dumps(CHECKS, indent=1))
    a = OUT / "figure_audit.csv"
    prev = pd.read_csv(a) if a.exists() else pd.DataFrame(columns=["figure", "files", "sources", "note"])
    new = pd.DataFrame(fm.AUDIT)
    pd.concat([prev[~prev.figure.isin(new.figure)], new]).sort_values("figure").to_csv(a, index=False)


if __name__ == "__main__":
    main()
