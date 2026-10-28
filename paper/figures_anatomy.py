"""Anatomical figures of the revised manuscript, drawn from the real mesh and the saved landmarks of
held-out vertebra S03 L4 (fold holding out S03; regions calibrated on S01 and S02; atlas S01 L1).

  fig_qualitative_insets   Fig. 7: expert reference, GCZ, and GCZ-RL with magnified insets of the
                           classes whose correction is active in this fold, and the per-class errors
                           of GCZ, GCZ-RL, and SLD on this vertebra
  fig_regions_margins      Fig. 8: expert attachment footprints and annotations, and GCZ-RL landmarks;
                           ALL and PLL with close-ups of the superior and inferior body margins.
                           Also Supplementary Fig. S8 (fig_search_regions): the calibrated GCZ search
                           regions, the masks used for projection, unchanged

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


# ------------------------------------------------------------------------------ verification table
def expert_part_labels(lig, pts):
    """Annotation sub-group (file suffix, e.g. CL_S_R) of every merged expert point of a class."""
    import re
    from scipy.spatial import cKDTree
    from data_io import LIGAMENT_SUFFIXES, landmark_dir, read_markup
    tags, coords = [], []
    for p in sorted(landmark_dir(*KEY).iterdir()):
        for suf in LIGAMENT_SUFFIXES[lig]:
            if re.match(rf".*_{re.escape(suf)}(_\d+)?\.mrk\.json$", p.name):
                q = read_markup(p)
                if q is not None:
                    coords.append(np.atleast_2d(q))
                    tags += [suf] * len(np.atleast_2d(q))
    d, idx = cKDTree(np.vstack(coords)).query(np.atleast_2d(pts))
    assert d.max() < 1e-6
    return [tags[i] for i in idx]


def visible(mesh, p, view, offset=0.4):
    """True if nothing of the mesh lies between the point and a camera of the given view."""
    vec = np.asarray(fs.VIEWS[view][0], float)
    start = np.asarray(p, float) + vec * offset
    hits, _ = mesh.ray_trace(start, start + vec * 400.0)
    return len(hits) == 0


def verify_landmarks():
    """Every landmark of S03 L4 drawn in Fig. 7: class, side, sub-group of its matched reference,
    coordinates of GCZ, GCZ-RL, and their Hungarian-matched references, errors, region membership,
    and visibility in the posterior and left lateral views. Recomputed class errors must equal the
    metrics file of the primary run."""
    import pyvista as pv
    from scipy.optimize import linear_sum_assignment
    from scipy.spatial import cKDTree
    v, f = fm.mesh(KEY)
    mesh = pv.PolyData(v, np.hstack([np.full((len(f), 1), 3), f]).ravel())
    gt = fm.cases.landmarks(*KEY)
    g0, rl = fm.preds("gcz")[KEY], fm.preds("gczrl")[KEY]
    check_saved("gcz", g0, fs.LIGAMENTS)
    check_saved("gczrl", rl, fs.LIGAMENTS)
    _, thr = fm.calibration()
    sel, proj = fm.region_masks(v, thr, third_fallback=True), fm.region_masks(v, thr, third_fallback=False)
    kd = cKDTree(v)
    met = pd.read_csv(fm.P2 / "per_vertebra_per_ligament_metrics.csv")
    met = met[met.case == f"{KEY[0]}_{KEY[1]}"].set_index(["method", "ligament"]).mae_mm
    xc = 0.5 * (v[:, 0].min() + v[:, 0].max())
    rows, vis_ref = [], []
    for lig in fs.LIGAMENTS:
        g = np.atleast_2d(gt[lig])
        part = expert_part_labels(lig, g)
        for k, p in enumerate(g):
            vis_ref.append({"ligament": lig, "part": part[k], "x": p[0], "y": p[1], "z": p[2],
                            "visible_posterior": visible(mesh, p, "posterior"), "visible_left": visible(mesh, p, "left")})
        match = {}
        for name, P in (("GCZ", g0), ("GCZ-RL", rl)):
            P = np.atleast_2d(P[lig])
            d = np.linalg.norm(P[:, None] - g[None], axis=2)
            ri, ci = linear_sum_assignment(d)
            match[name] = dict(zip(ri, ci))
            err = float(d[ri, ci].mean())
            assert abs(err - met[(name, lig)]) < 1e-9, (name, lig, err, met[(name, lig)])
        A, B = np.atleast_2d(g0[lig]), np.atleast_2d(rl[lig])
        for i in range(len(A)):
            ra, rb = g[match["GCZ"][i]], g[match["GCZ-RL"][i]]
            rows.append({
                "ligament": lig, "point_index": i, "side": "left" if B[i, 0] > xc else "right",
                **{f"gcz_{c}": A[i, j] for j, c in enumerate("xyz")}, **{f"gczrl_{c}": B[i, j] for j, c in enumerate("xyz")},
                "correction_mm": float(np.linalg.norm(B[i] - A[i])),
                "gcz_ref_part": part[match["GCZ"][i]], **{f"gcz_ref_{c}": ra[j] for j, c in enumerate("xyz")},
                "gcz_error_mm": float(np.linalg.norm(A[i] - ra)),
                "gczrl_ref_part": part[match["GCZ-RL"][i]], **{f"gczrl_ref_{c}": rb[j] for j, c in enumerate("xyz")},
                "gczrl_error_mm": float(np.linalg.norm(B[i] - rb)),
                "gcz_dist_to_mesh_mm": float(kd.query(A[i])[0]), "gczrl_dist_to_mesh_mm": float(kd.query(B[i])[0]),
                "gcz_in_selection_region": bool(sel[lig][kd.query(A[i])[1]]),
                "gcz_in_projection_region": bool(proj[lig][kd.query(A[i])[1]]),
                "gczrl_in_projection_region": bool(proj[lig][kd.query(B[i])[1]]),
                "gcz_visible_posterior": visible(mesh, A[i], "posterior"), "gczrl_visible_posterior": visible(mesh, B[i], "posterior"),
                "gcz_visible_left": visible(mesh, A[i], "left"), "gczrl_visible_left": visible(mesh, B[i], "left")})
    tab = pd.DataFrame(rows)
    tab.round(3).to_csv(OUT / "fig7_landmark_verification_S03_L4.csv", index=False)
    pd.DataFrame(vis_ref).round(3).to_csv(OUT / "fig7_reference_visibility_S03_L4.csv", index=False)
    CHECKS["fig7_class_errors_recomputed_equal_metrics_file"] = True
    CHECKS["fig7_all_predictions_on_mesh_vertices"] = bool((tab[["gcz_dist_to_mesh_mm", "gczrl_dist_to_mesh_mm"]] < 1e-4).all().all())
    return tab, pd.DataFrame(vis_ref)


# ------------------------------------------------------------------------------ Fig. 7
def marker_visible(mesh, p, view, r):
    """A marker of radius r is visible if the ray towards the camera is free from its centre or from
    one of four points at 0.7 r around it (perpendicular to the view direction)."""
    vec = np.asarray(fs.VIEWS[view][0], float)
    a = np.cross(vec, [0.0, 0.0, 1.0])
    a = a / np.linalg.norm(a) if np.linalg.norm(a) > 1e-6 else np.array([1.0, 0.0, 0.0])
    b = np.cross(vec, a)
    for o in (np.zeros(3), a, -a, b, -b):
        if visible(mesh, np.asarray(p, float) + 0.7 * r * o, view):
            return True
    return False


def fig_qualitative_insets():
    """Fig. 7: one posterior camera. The overview (a) and every close-up (b-e) are the same render,
    so each close-up is exactly the pixels inside its frame on the overview."""
    import pyvista as pv
    v, f = fm.mesh(KEY)
    mesh = pv.PolyData(v, np.hstack([np.full((len(f), 1), 3), f]).ravel())
    gt = fm.cases.landmarks(*KEY)
    g0, rl = fm.preds("gcz")[KEY], fm.preds("gczrl")[KEY]
    tab, _ = verify_landmarks()
    met = pd.read_csv(fm.P2 / "per_vertebra_per_ligament_metrics.csv")
    met = met[met.case == f"{KEY[0]}_{KEY[1]}"].pivot(index="ligament", columns="method", values="mae_mm")
    xc = 0.5 * (v[:, 0].min() + v[:, 0].max())
    ligs = ["ITL", "CL", "ISL", "SSL"]                      # the classes whose correction is active in this fold
    pairs = [p for p in fm.correction_pairs(g0, rl, ligs) if np.linalg.norm(np.subtract(*p)) > 0.3]
    lays = [layer("spheres", stack(gt, ligs), EXPERT_C, R_EXPERT), layer("cubes", stack(g0, ligs), GCZ_C, E_GCZ),
            {"kind": "lines", "pairs": pairs, "color": SEG_C, "radius": 0.15},
            layer("spheres", stack(rl, ligs), RL_C, R_RL)]
    A = lambda d, l: np.atleast_2d(d[l])  # noqa: E731
    side = lambda P, s: P[(P[:, 0] > xc) == (s == "left")]  # noqa: E731  (LPS: +x is the patient's left)
    # right CL around the superior facet joint: the CL_S_R references and the landmarks paired with them
    cl = tab[(tab.ligament == "CL") & (tab.gczrl_ref_part == "CL_S_R")]
    parts = np.array(expert_part_labels("CL", A(gt, "CL")))
    cl_pts = np.vstack([A(gt, "CL")[parts == "CL_S_R"], cl[["gcz_x", "gcz_y", "gcz_z"]].to_numpy(),
                        cl[["gczrl_x", "gczrl_y", "gczrl_z"]].to_numpy()])
    sp_pts = np.vstack([A(gt, "ISL"), A(gt, "SSL"), A(g0, "ISL"), A(g0, "SSL"), A(rl, "ISL"), A(rl, "SSL")])
    itl = {s: np.vstack([side(A(gt, "ITL"), s), side(A(g0, "ITL"), s), side(A(rl, "ITL"), s)]) for s in ("right", "left")}
    sets = [cl_pts, sp_pts, itl["right"], itl["left"]]
    img, xy, pxmm = render(v, f, "posterior", lays, size=3200, project=np.vstack(sets))
    cut = np.cumsum([0] + [len(s) for s in sets])
    P = [xy[cut[i]:cut[i + 1]] for i in range(len(sets))]
    boxes = {"b": box_around(P[0], pxmm, 3.0, aspect=1.0), "c": box_around(P[1], pxmm, 2.5)}
    # ITL: the whole transverse process, from 5 mm beyond the tip to 30 mm towards the midline
    mid_px = render(v, f, "posterior", (), size=3200, project=[[xc, v[:, 1].mean(), v[:, 2].mean()]])[1][0, 0]
    for letter, s, pts in (("d", "right", P[2]), ("e", "left", P[3])):
        tip = pts[0]                                         # the expert point is the first row of the set
        inward = np.sign(mid_px - tip[0])
        x_out, x_in = tip[0] - inward * 5 * pxmm, tip[0] + inward * 30 * pxmm
        yc = 0.5 * (pts[:, 1].min() + pts[:, 1].max())
        boxes[letter] = np.array([min(x_out, x_in), yc - 10 * pxmm, max(x_out, x_in), yc + 10 * pxmm])
        assert (pts[:, 0] >= boxes[letter][0]).all() and (pts[:, 0] <= boxes[letter][2]).all()
    for k, s in zip("bcde", sets):
        q = P["bcde".index(k)]
        assert ((q >= boxes[k][:2]) & (q <= boxes[k][2:])).all(), k     # every defining point is inside its frame
    CHECKS["fig7_closeups_are_crops_of_the_overview"] = True
    CHECKS["fig7_camera"] = "posterior (+y towards the camera, +z up); the patient's right is on the right"
    CHECKS["fig7_reference_markers_visible_in_closeups"] = {
        "b (CL_S_R)": f"{sum(marker_visible(mesh, p, 'posterior', R_EXPERT) for p in A(gt, 'CL')[parts == 'CL_S_R'])}"
                      f"/{int((parts == 'CL_S_R').sum())}",
        "c (ISL, SSL)": f"{sum(marker_visible(mesh, p, 'posterior', R_EXPERT) for p in np.vstack([A(gt, 'ISL'), A(gt, 'SSL')]))}"
                        f"/{len(A(gt, 'ISL')) + len(A(gt, 'SSL'))}",
        "d, e (ITL)": f"{sum(marker_visible(mesh, p, 'posterior', R_EXPERT) for p in A(gt, 'ITL'))}/{len(A(gt, 'ITL'))}"}
    CHECKS["fig7_cl_closeup_landmarks"] = cl[["point_index", "side", "gcz_error_mm", "gczrl_error_mm"]].round(3).to_dict("records")

    fig = plt.figure(figsize=(7.2, 5.0))
    outer = fig.add_gridspec(2, 1, height_ratios=[1.15, 1.0], hspace=0.22)
    fig.subplots_adjust(left=0.01, right=0.99, top=0.91, bottom=0.075)
    top = outer[0].subgridspec(1, 3, width_ratios=[2.15, 1, 1], wspace=0.06)
    cb = union(content_box(img, pad=24), *[b + [-90, -110, 20, 20] for b in boxes.values()])
    ax = fig.add_subplot(top[0, 0])
    show(ax, crop(img, cb), "(a) S03 L4, posterior view")
    for k, b in boxes.items():
        frame_outside(ax, b, cb[:2], k, others=[o for j, o in boxes.items() if j != k])
    show(fig.add_subplot(top[0, 1]), crop(img, boxes["b"]), "(b) CL, right facet joint", boxed=True)
    show(fig.add_subplot(top[0, 2]), crop(img, boxes["c"]), "(c) ISL and SSL", boxed=True)
    bot = outer[1].subgridspec(1, 3, width_ratios=[1.35, 1.35, 1.6], wspace=0.1)
    show(fig.add_subplot(bot[0, 0]), crop(img, boxes["d"]), "(d) ITL, right transverse process", boxed=True)
    show(fig.add_subplot(bot[0, 1]), crop(img, boxes["e"]), "(e) ITL, left transverse process", boxed=True)

    ax = fig.add_subplot(bot[0, 2])
    pos = ax.get_position()
    ax.set_position([pos.x0 + 0.05, pos.y0 + 0.02, pos.width - 0.05, pos.height - 0.02])
    order = fs.LIGAMENTS[::-1]
    for i, l in enumerate(order):
        a0, a1, s = met.loc[l, "GCZ"], met.loc[l, "GCZ-RL"], met.loc[l, "SLD"]
        if abs(a0 - a1) > 1e-9:
            ax.annotate("", xy=(a1, i), xytext=(a0, i),
                        arrowprops=dict(arrowstyle="-|>", color=fs.NEUTRAL, lw=0.9, shrinkA=3.5, shrinkB=3.5,
                                        mutation_scale=7))
        ax.plot(a0, i, "s", ms=4.4, color=GCZ_C, zorder=3)
        ax.plot(s, i, "^", ms=4.6, color=SLD_C, zorder=4)
        ax.plot(a1, i, "o", ms=4.6, color=RL_C, zorder=5)
    ax.set_yticks(range(len(order)), order)
    ax.set_ylim(-0.6, len(order) - 0.4)
    ax.set_xlim(0, 8.6)
    ax.set_xlabel("Class MAE on S03 L4 (mm)")
    ax.xaxis.grid(True, color=fs.GRID, linewidth=0.5)
    ax.set_axisbelow(True)
    ax.tick_params(axis="y", length=0)
    ax.set_title("(f) Error by class", loc="left", fontsize=7.5, fontweight="bold", pad=3)
    mean = met.mean()
    CHECKS["fig7_macro_mae_S03_L4"] = {k: round(float(mean[k]), 3) for k in ("GCZ", "GCZ-RL", "SLD")}
    CHECKS["fig7_class_mae_S03_L4"] = {l: {k: round(float(met.loc[l, k]), 3) for k in ("GCZ", "GCZ-RL", "SLD")}
                                       for l in fs.LIGAMENTS}
    CHECKS["fig7_itl_points"] = {r.side: {"gcz_error_mm": round(r.gcz_error_mm, 3), "gczrl_error_mm": round(r.gczrl_error_mm, 3),
                                          "gcz_in_projection_region": bool(r.gcz_in_projection_region),
                                          "gcz_in_selection_region": bool(r.gcz_in_selection_region)}
                                 for r in tab[tab.ligament == "ITL"].itertuples()}
    from matplotlib.lines import Line2D
    from matplotlib.patches import Rectangle
    h = [Line2D([], [], ls="", marker="o", color=EXPERT_C, ms=3.4, label="Expert reference"),
         Line2D([], [], ls="", marker="s", color=GCZ_C, ms=4.6, label="GCZ (atlas transfer)"),
         Line2D([], [], ls="", marker="o", color=RL_C, ms=4.8, label="GCZ-RL (corrected)"),
         Line2D([], [], color=SEG_C, lw=0.9, label="Applied correction"),
         Line2D([], [], ls="", marker="^", color=SLD_C, ms=4.8, label="SLD (panel f only)"),
         Rectangle((0, 0), 1, 1, fill=False, ec=fs.INK, lw=0.8, label="Close-up frame")]
    fig.legend(handles=h, loc="upper center", ncol=6, frameon=False, bbox_to_anchor=(0.5, 0.998), fontsize=6.7,
               handletextpad=0.4, columnspacing=1.1)
    files = fs.save(fig, OUT, "fig_qualitative_insets", vector=False)
    fm.audit("Fig. 7 (new)", files, [fm.P2 / "predictions_gcz.csv", fm.P2 / "predictions_gczrl.csv",
                                     fm.P2 / "per_vertebra_per_ligament_metrics.csv", "raw mesh and annotations of S03 L4",
                                     "fig7_landmark_verification_S03_L4.csv"])
    return files


def frame_outside(ax, b, origin, label, color=fs.INK, others=()):
    """Thin unfilled close-up frame with its letter just outside the first of its four corners that
    does not fall inside another frame."""
    from matplotlib.patches import Rectangle
    x0, y0, x1, y1 = b
    ox, oy = origin
    ax.add_patch(Rectangle((x0 - ox, y0 - oy), x1 - x0, y1 - y0, fill=False, ec=color, lw=0.8))
    d = 40
    corners = [((x0 - d, y0 - d), "right", "bottom"), ((x1 + d, y0 - d), "left", "bottom"),
               ((x0 - d, y1 + d), "right", "top"), ((x1 + d, y1 + d), "left", "top")]
    inside = lambda q: any(o[0] <= q[0] <= o[2] and o[1] <= q[1] <= o[3] for o in others)  # noqa: E731
    (qx, qy), ha, va = next((c for c in corners if not inside(c[0])), corners[0])
    ax.text(qx - ox + (d - 8) * (1 if ha == "right" else -1), qy - oy + (d - 8) * (1 if va == "bottom" else -1),
            label, ha=ha, va=va, fontsize=7.4, fontweight="bold", color=color)


FOOT_MM = 2.5                                         # footprint: surface within 2.5 mm of an expert point
FOOT_S, FOOT_I, FOOT = "#86c3e6", "#c3b1e2", "#9ed4c6"  # superior margin, inferior margin, other classes


def fig_regions_margins():
    from scipy.spatial import cKDTree
    v, f = fm.mesh(KEY)
    gt = fm.cases.landmarks(*KEY)
    rl = fm.preds("gczrl")[KEY]
    check_saved("gczrl", rl, fs.LIGAMENTS)
    _, thr = fm.calibration()
    masks = fm.region_masks(v, thr)
    ref = mv.region_mask(v, "ALL", mv.thresholds("Verse581"))
    CHECKS["fig8_masks_equal_margin_verification"] = bool(np.array_equal(ref, masks["ALL"]))
    walls = mv.body_walls(v)
    kd = cKDTree(v)
    rgb_of = lambda c: (np.array(fs.matplotlib.colors.to_rgb(c)) * 255).astype(np.uint8)  # noqa: E731
    bone = rgb_of(fs.BONE)

    def footprint(pts):
        m = np.zeros(len(v), bool)
        for idx in kd.query_ball_point(np.atleast_2d(pts), FOOT_MM):
            m[idx] = True
        return m

    def paint(parts):
        rgb = np.tile(bone, (len(v), 1))
        for m, c in parts:
            rgb[m] = rgb_of(c)
        return rgb

    foot = {}
    for lig in ("ALL", "PLL"):
        pp = mv.parts(*KEY, lig)
        foot[lig] = {"S": footprint(pp["S"]), "I": footprint(pp["I"])}
    for lig in ("LF", "ISL", "SSL", "ITL", "CL"):
        foot[lig] = footprint(gt[lig])
    # every footprint vertex lies within FOOT_MM of an expert point of its class; footprints are not masks
    CHECKS["fig8_footprint_vertices_within_mm_of_expert"] = FOOT_MM
    CHECKS["fig8_footprint_share_of_mask_area_check"] = {
        l: round(float(((foot[l]["S"] | foot[l]["I"]) if l in ("ALL", "PLL") else foot[l]).sum() / max(masks[l].sum(), 1)), 3)
        for l in fs.LIGAMENTS}

    def letter(ax, s):
        ax.text(0.025, 0.97, s, transform=ax.transAxes, ha="left", va="top", fontsize=7.5, fontweight="bold",
                bbox=dict(fc="white", ec="none", pad=0.6, alpha=0.85))

    r_e, r_p = 0.7, 1.0
    exp_l = lambda l, r=r_e: layer("spheres", gt[l], EXPERT_C, r)  # noqa: E731
    rl_l = lambda l: layer("spheres", rl[l], RL_C, r_p)  # noqa: E731

    fig = plt.figure(figsize=(7.2, 5.75))
    outer = fig.add_gridspec(3, 1, height_ratios=[1.2, 1.45, 1.45], hspace=0.3)
    fig.subplots_adjust(left=0.075, right=0.995, top=0.915, bottom=0.005)

    # ---- Fig. 8: expert footprints, expert points, GCZ-RL landmarks
    STRIP_W, STRIP_H = 46.0, 21.0                 # close-up size, mm
    fb = mv.body_only(v, f, walls["y_post"] + 0.5)
    for r, (lig, view, faces) in enumerate((("ALL", "anterior", f), ("PLL", "posterior", fb))):
        row = outer[r].subgridspec(1, 3, width_ratios=[1.0, 1.25, 1.25], wspace=0.05)
        pp = mv.parts(*KEY, lig)
        body = v[np.unique(faces)]
        center = 0.5 * (body.min(0) + body.max(0))
        scale = 0.5 * max(np.ptp(body[:, 0]), np.ptp(body[:, 2])) * 1.1
        rgb = paint([(foot[lig]["S"], FOOT_S), (foot[lig]["I"], FOOT_I)])
        img, xy, pxmm = render(v, faces, view, [exp_l(lig), rl_l(lig)], rgb=rgb, center=center, scale=scale,
                               project=np.vstack([pp["S"].mean(0), pp["I"].mean(0)]))
        boxes = [np.r_[c - [STRIP_W / 2 * pxmm, STRIP_H / 2 * pxmm], c + [STRIP_W / 2 * pxmm, STRIP_H / 2 * pxmm]]
                 for c in xy]
        cbx = union(content_box(img, pad=15), *[b + [-8, -8, 8, 8] for b in boxes])
        ax = fig.add_subplot(row[0, 0])
        show(ax, crop(img, cbx))
        letter(ax, "ad"[r])
        for b in boxes:
            frame(ax, b, cbx[:2])
        ax.text(-0.04, 0.5, f"{lig}\n({'anterior' if lig == 'ALL' else 'posterior'})", transform=ax.transAxes,
                rotation=90, ha="right", va="center", fontsize=7.2, fontweight="bold")
        for c, b in enumerate(boxes):
            a = fig.add_subplot(row[0, c + 1])
            show(a, crop(img, b), boxed=True)
            letter(a, ("bc" if lig == "ALL" else "ef")[c])
        if r == 0:
            for a, t in zip(fig.axes[-3:], ["Whole view with close-up frames", "Superior margin (close-up)",
                                            "Inferior margin (close-up)"]):
                a.set_title(t, loc="left", fontsize=7.3, fontweight="bold", pad=3)

    xc = 0.5 * (v[:, 0].min() + v[:, 0].max())
    left = lambda p: np.atleast_2d(p)[np.atleast_2d(p)[:, 0] > xc]  # noqa: E731  (LPS: +x is the left)
    spec = {"LF": ("superior", np.atleast_2d(gt["LF"]), "LF, superior view", None),
            "ISL": ("left", np.vstack([gt["ISL"], gt["SSL"]]), "ISL, left lateral", None),
            "SSL": ("posterior", np.atleast_2d(gt["SSL"]), "SSL, tip (posterior)", 24),
            "ITL": ("left", left(gt["ITL"]), "ITL, left tip (lateral)", 24),
            "CL": ("posterior", left(gt["CL"]), "CL, left side (posterior)", None)}
    row = outer[2].subgridspec(1, 5, wspace=0.05)
    for c, lig in enumerate(["LF", "ISL", "SSL", "ITL", "CL"]):
        view, focus, title, mn = spec[lig]
        img, xy, pxmm = render(v, f, view, [exp_l(lig, 0.75 if lig in ("SSL", "ITL") else r_e), rl_l(lig)],
                               rgb=paint([(foot[lig], FOOT)]), project=focus)
        box = box_around(xy, pxmm, 3.0, aspect=1.0, min_mm=mn)
        ax = fig.add_subplot(row[0, c])
        show(ax, crop(img, box), boxed=True)
        letter(ax, "ghijk"[c])
        ax.set_title(title, loc="left", fontsize=7.0, fontweight="bold", pad=3)
        spec[lig] = (view, focus, title, mn, box)
        if c == 0:
            ax.text(-0.06, 0.5, "Expert footprint,\nexpert points,\nGCZ-RL", transform=ax.transAxes, rotation=90,
                    ha="right", va="center", fontsize=6.8)
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch, Rectangle
    h = [Line2D([], [], ls="", marker="o", color=EXPERT_C, ms=3.4, label="Expert annotation"),
         Line2D([], [], ls="", marker="o", color=RL_C, ms=4.8, label="GCZ-RL landmark"),
         Patch(fc=FOOT_S, ec="none", label="Footprint, superior margin"),
         Patch(fc=FOOT_I, ec="none", label="Footprint, inferior margin"),
         Patch(fc=FOOT, ec="none", label="Footprint, other classes"),
         Rectangle((0, 0), 1, 1, fill=False, ec=fs.INK, lw=0.8, label="Close-up frame")]
    fig.legend(handles=h, loc="upper center", ncol=3, frameon=False, bbox_to_anchor=(0.5, 1.0), fontsize=6.7,
               handletextpad=0.4, columnspacing=1.2)
    CHECKS["fig8_rl_inside_region"] = {l: bool(np.isin(np.round(np.atleast_2d(rl[l]), 4),
                                                       np.round(v[masks[l]], 4)).all(axis=1).all())
                                       for l in fs.LIGAMENTS}
    CHECKS["fig8_closeups_are_crops_of_the_whole_view"] = True
    files = fs.save(fig, OUT, "fig_regions_margins", vector=False)
    fm.audit("Fig. 8 (new)", files, [fm.P2 / "predictions_gczrl.csv", "raw mesh and annotations of S03 L4"])

    # ---- Supplementary Fig. S8: the calibrated GCZ search regions, exactly as computed, with the expert points
    fig = plt.figure(figsize=(7.2, 1.68))
    fig.subplots_adjust(left=0.075, right=0.995, top=0.86, bottom=0.01)
    row = fig.add_gridspec(1, 7, wspace=0.05)
    grey = rgb_of(REGION_C)
    for c, lig in enumerate(fs.LIGAMENTS):
        rgb = np.tile(bone, (len(v), 1))
        rgb[masks[lig]] = grey
        if lig in ("ALL", "PLL"):
            faces = f if lig == "ALL" else fb
            body = v[np.unique(faces)]
            img, _, _ = render(v, faces, "anterior" if lig == "ALL" else "posterior", [exp_l(lig, 0.6)], rgb=rgb,
                               center=0.5 * (body.min(0) + body.max(0)),
                               scale=0.5 * max(np.ptp(body[:, 0]), np.ptp(body[:, 2])) * 1.1)
            box = content_box(img, pad=10)
            box = box_around(np.array([[box[0], box[1]], [box[2], box[3]]]), 1.0, 0.0, aspect=1.0)
        else:
            view, focus, _, mn, box = spec[lig]
            img, _, _ = render(v, f, view, [exp_l(lig, 0.75 if lig in ("SSL", "ITL") else 0.6)], rgb=rgb)
        ax = fig.add_subplot(row[0, c])
        show(ax, crop(img, box), boxed=True)
        letter(ax, "abcdefg"[c])
        ax.set_title(lig, loc="left", fontsize=7.0, fontweight="bold", pad=3)
        if c == 0:
            ax.text(-0.06, 0.5, "Calibrated GCZ\nsearch region\n(actual mask)", transform=ax.transAxes,
                    rotation=90, ha="right", va="center", fontsize=6.6)
        # the gray vertices of the render input are exactly the mask vertices
        assert np.array_equal(np.all(rgb == grey, axis=1), masks[lig])
    h = [Line2D([], [], ls="", marker="o", color=EXPERT_C, ms=3.4, label="Expert annotation"),
         Patch(fc=REGION_C, ec="none", label="Calibrated GCZ search region (the actual mask used for projection, unchanged)")]
    fig.legend(handles=h, loc="upper center", ncol=2, frameon=False, bbox_to_anchor=(0.5, 1.0), fontsize=6.7,
               handletextpad=0.4, columnspacing=1.4)
    CHECKS["figS8_gray_vertices_equal_masks"] = True
    s8 = fs.save(fig, OUT, "fig_search_regions", vector=False)
    fm.audit("Fig. S8 (new)", s8, ["gcz.build_thresholds / candidate_vertices (fold holding out S03)",
                                   "raw mesh and annotations of S03 L4"])
    return files + s8


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
