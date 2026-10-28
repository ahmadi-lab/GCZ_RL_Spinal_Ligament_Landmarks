"""Anatomical verification of the GCZ candidate regions against the expert attachment annotations.

For every held-out vertebra of the primary evaluation, the regions are rebuilt exactly as in the
method from the training subjects of its fold (gcz.build_thresholds, gcz.candidate_vertices) and
compared with the expert points of the same class. Nothing is changed: the masks, predictions,
and metrics are those of the evaluated pipeline.

Per vertebra and class:
  band_*            calibrated bands (lo-hi) of y_n, |x_n-0.5|, z_n
  inside_*_pct      expert points whose nearest vertex lies in the selection / projection region
  region_area_pct   share of the vertebral surface area in the projection region
  footprint_area_pct  share of the surface within FOOT_MM of an expert point of the class
  region_on_footprint_pct  share of the region area within FOOT_MM of an expert point
  footprint_in_region_pct  share of the footprint area inside the region
  overlap_pct       share of the region area that also belongs to another class region
  extent_*          extent of the region and of the expert points along x (lateral), y (AP), z (SI), mm
  rl_inside_pct     GCZ-RL landmarks that are vertices of the projection region

Also checks that the rendering colours exactly the mask vertices, and draws the figure that shows
expert footprints, candidate regions, and GCZ-RL landmarks separately (held-out S03 L4).

Usage: python paper/region_verification.py
Outputs: paper/outputs/region_verification/*.csv, summary.json, fig_regions_footprints.{pdf,png}
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

OUT = REPO / "paper" / "outputs" / "region_verification"
P2 = REPO / "part2_patient_loso_14" / "outputs"
COHORT = {"Verse005": ["L1", "L2", "L3", "L4", "L5"], "Verse561": ["L1", "L2", "L3", "L4"],
          "Verse581": ["L1", "L2", "L3", "L4", "L5"]}
FOOT_MM = 2.5
LIGS = fs.LIGAMENTS
VIEW_OF = {"ALL": "anterior", "PLL": "superior", "LF": "superior", "ISL": "left", "SSL": "left",
           "ITL": "posterior", "CL": "posterior"}
fs.VIEWS.setdefault("anterior", ((0, -1, 0), (0, 0, 1)))
fs.VIEW_TITLE.setdefault("anterior", "Anterior view")


def vertex_area(v, f):
    tri = v[f]
    a = 0.5 * np.linalg.norm(np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0]), axis=1)
    va = np.zeros(len(v))
    for c in range(3):
        np.add.at(va, f[:, c], a / 3.0)
    return va


def mask_of(v, zone, kd):
    d, idx = kd.query(np.asarray(zone, float))
    m = np.zeros(len(v), bool)
    m[idx[d < 1e-6]] = True
    return m


def fold_thresholds(held):
    train = {(s, l): {"verts": cases.mesh(s, l, "vtk"), "gt": cases.landmarks(s, l)}
             for s in COHORT if s != held for l in COHORT[s]}
    return gcz.build_thresholds(train)


def per_vertebra():
    rl_pred = table_from_csv(P2 / "predictions_gczrl.csv")
    rows, masks_s03l4 = [], None
    for held in COHORT:
        thr = fold_thresholds(held)
        for lv in COHORT[held]:
            v, f = read_vtk_polydata(RAW_ROOT / held / f"{held}_{lv}.vtk")
            v = np.asarray(v, float)
            va = vertex_area(v, f)
            kd = cKDTree(v)
            gt = cases.landmarks(held, lv)
            proj = {l: mask_of(v, gcz.candidate_vertices(v, l, thr, third_fallback=False), kd) for l in LIGS}
            sel = {l: mask_of(v, gcz.candidate_vertices(v, l, thr, third_fallback=True), kd) for l in LIGS}
            if (held, lv) == ("Verse581", "L4"):
                masks_s03l4 = (thr, proj)
            for l in LIGS:
                g = gt.get(l)
                if g is None or len(g) == 0:
                    continue
                g = np.atleast_2d(g)
                nn = kd.query(g)[1]
                foot = np.zeros(len(v), bool)
                for idx in kd.query_ball_point(g, FOOT_MM):
                    foot[idx] = True
                reg = proj[l]
                others = np.zeros(len(v), bool)
                for o in LIGS:
                    if o != l:
                        others |= proj[o]
                p = rl_pred.get((held, lv), {}).get(l)
                p_inside = None
                if p is not None and len(p):
                    p_inside = float(100 * (cKDTree(v[reg]).query(np.atleast_2d(p))[0] < 1e-4).mean()) if reg.any() else 0.0
                ext = lambda pts: np.ptp(pts, axis=0) if len(pts) > 1 else np.zeros(3)  # noqa: E731
                er, eg = ext(v[reg]), ext(g)
                rows.append({
                    "held_out": SUBJECT_CODE[held], "case": f"{SUBJECT_CODE[held]} {lv}", "ligament": l, "n_ref": len(g),
                    **{f"band_{q}": f"{thr[l][q]['lo']:.2f}-{thr[l][q]['hi']:.2f}" for q in ("y_n", "x_half", "z_n")},
                    "inside_selection_pct": 100 * sel[l][nn].mean(), "inside_projection_pct": 100 * reg[nn].mean(),
                    "region_vertices": int(reg.sum()), "region_area_pct": 100 * va[reg].sum() / va.sum(),
                    "footprint_area_pct": 100 * va[foot].sum() / va.sum(),
                    "region_on_footprint_pct": 100 * va[reg & foot].sum() / max(va[reg].sum(), 1e-12),
                    "footprint_in_region_pct": 100 * va[reg & foot].sum() / max(va[foot].sum(), 1e-12),
                    "overlap_pct": 100 * va[reg & others].sum() / max(va[reg].sum(), 1e-12),
                    "extent_region_x": er[0], "extent_region_y": er[1], "extent_region_z": er[2],
                    "extent_expert_x": eg[0], "extent_expert_y": eg[1], "extent_expert_z": eg[2],
                    "rl_inside_pct": p_inside})
    return pd.DataFrame(rows), masks_s03l4


def summarize(df):
    agg = []
    for l in LIGS:
        x = df[df.ligament == l]
        agg.append({"ligament": l, "vertebrae": len(x), "n_ref": int(x.n_ref.sum()),
                    "inside_selection_pct": round(100 * (x.inside_selection_pct * x.n_ref).sum() / 100 / x.n_ref.sum(), 1),
                    "region_area_pct": round(x.region_area_pct.median(), 1),
                    "footprint_area_pct": round(x.footprint_area_pct.median(), 1),
                    "region_to_footprint_ratio": round((x.region_area_pct / x.footprint_area_pct).median(), 1),
                    "region_on_footprint_pct": round(x.region_on_footprint_pct.median(), 1),
                    "footprint_in_region_pct": round(x.footprint_in_region_pct.median(), 1),
                    "overlap_pct": round(x.overlap_pct.median(), 1),
                    "extent_region_xyz_mm": "/".join(f"{x[f'extent_region_{a}'].median():.0f}" for a in "xyz"),
                    "extent_expert_xyz_mm": "/".join(f"{x[f'extent_expert_{a}'].median():.0f}" for a in "xyz"),
                    "rl_inside_pct": round(x.rl_inside_pct.mean(), 1)})
    return pd.DataFrame(agg)


def render_check(masks):
    """The rendering colours exactly the vertices of the region masks."""
    v, f = read_vtk_polydata(RAW_ROOT / "Verse581" / "Verse581_L4.vtk")
    _, proj = masks
    ok = True
    for l in LIGS:
        rgb = fs.region_colors(np.asarray(v, float), {l: proj[l]}, order=(l,), mix=0.6)
        bone = (np.array(fs.matplotlib.colors.to_rgb(fs.BONE)) * 255).astype(np.uint8)
        coloured = ~np.all(rgb == bone, axis=1)
        ok &= bool(np.array_equal(coloured, proj[l]))
    return ok


def figure(masks):
    fs.setup()
    plt = fs.plt
    key = ("Verse581", "L4")
    v, f = read_vtk_polydata(RAW_ROOT / key[0] / f"{key[0]}_{key[1]}.vtk")
    v = np.asarray(v, float)
    gt = cases.landmarks(*key)
    rl = table_from_csv(P2 / "predictions_gczrl.csv")[key]
    _, proj = masks
    rows = ["Expert attachment\nannotations", "GCZ candidate\nregion", "GCZ-RL\nlandmarks"]
    fig, axes = plt.subplots(3, len(LIGS), figsize=(7.2, 3.6))
    for c, l in enumerate(LIGS):
        view = VIEW_OF[l]
        g = np.atleast_2d(gt[l])
        # very small regions (SSL, ITL: tens of vertices) are also drawn as markers on the same vertices
        small = [{"kind": "spheres", "points": v[proj[l]], "color": fs.CLASS[l], "size": 12}] if proj[l].sum() < 300 else []
        imgs = [fs.render(v, f, view, [{"kind": "spheres", "points": g, "color": fs.CLASS[l], "size": 20}], zoom=1.25),
                fs.render(v, f, view, small, zoom=1.25, vertex_rgb=fs.region_colors(v, {l: proj[l]}, order=(l,), mix=0.6)),
                fs.render(v, f, view, [{"kind": "spheres", "points": g, "color": "#9a9a9a", "size": 9},
                                       {"kind": "spheres", "points": rl[l], "color": fs.CLASS[l], "size": 30}],
                          zoom=1.25)]
        for r_, img in enumerate(imgs):
            ax = axes[r_, c]
            fs.show(ax, img)
            if r_ == 0:
                ax.set_title(f"{l}\n({fs.VIEW_TITLE[view].replace(' view', '').lower()})", fontsize=7, pad=2)
            if c == 0:
                ax.set_ylabel(rows[r_], fontsize=7)
    fig.subplots_adjust(left=0.07, right=0.995, top=0.88, bottom=0.01, wspace=0.03, hspace=0.06)
    return fs.save(fig, OUT, "fig_regions_footprints", vector=False)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    df, masks = per_vertebra()
    df.round(2).to_csv(OUT / "region_verification_per_vertebra.csv", index=False)
    summ = summarize(df)
    summ.to_csv(OUT / "region_verification_by_class.csv", index=False)
    bands = df[df.held_out == "S03"].drop_duplicates("ligament")[["ligament", "band_y_n", "band_x_half", "band_z_n"]]
    bands.to_csv(OUT / "bands_fold_S03.csv", index=False)
    band = bands.set_index("ligament")
    rows = [f"{r.ligament} & {band.loc[r.ligament, 'band_y_n']} & {band.loc[r.ligament, 'band_x_half']} & "
            f"{band.loc[r.ligament, 'band_z_n']} & {r.inside_selection_pct:.1f} & {r.region_area_pct:.1f} & "
            f"{r.footprint_area_pct:.1f} & {r.region_to_footprint_ratio:.1f} & {r.footprint_in_region_pct:.1f} & "
            f"{r.region_on_footprint_pct:.1f} & {r.overlap_pct:.1f} \\\\".replace("-", "--")
            for r in summ.itertuples()]
    (OUT / "tab_region_verification.tex").write_text("\n".join(rows) + "\n", encoding="utf-8")
    render_ok = render_check(masks)
    files = figure(masks)
    (OUT / "summary.json").write_text(json.dumps({"render_colours_equal_masks": render_ok,
                                                  "footprint_buffer_mm": FOOT_MM,
                                                  "figure": [p.name for p in files]}, indent=2))
    print(summ.to_string(index=False))
    print(bands.to_string(index=False))
    print("rendering colours exactly the mask vertices:", render_ok)


if __name__ == "__main__":
    main()
