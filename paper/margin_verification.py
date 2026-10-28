"""Verification of the ALL and PLL candidate regions against the superior and inferior body margins.

For the 14 held-out vertebrae of the primary evaluation (regions rebuilt from the training
subjects of each fold, exactly as in the method) the script checks:
  1. orientation: superior (_S) expert points lie above inferior (_I) points along +z (LPS frame),
     by at least the body height;
  2. normalization: z_n of the code equals (z - z_min) / (z_max - z_min) over the whole vertebra;
  3. calibration: the z_n positions of the superior and inferior parts in each training vertebra,
     and the pooled band mu +/- 2 sigma that the method calibrates from them;
  4. masks: the share of the region surface, of the expert points, and of the GCZ and GCZ-RL
     landmarks that lies in the middle half between the body margins (s = 0.25-0.75, where s = 0
     at the inferior and 1 at the superior margin of the same wall, each margin taken as the mean
     height of the expert points of that part);
  5. what separate superior and inferior bands would change: the share of the region and of the
     saved landmarks inside the union of the two separate bands (a counterfactual check only).
In the figure (S03 L4) the margins are also measured on the mesh alone: the anterior wall (ALL) is
the front 5 mm of the body within 8 mm of the midline; the posterior wall (PLL) is the 3 mm in
front of the spinal canal, found as the first gap of at least 4 mm along y in the midline slab.
This simple detector assumes a vertebra that is not tilted in the sagittal plane, so it is used
for the figure case only, where it agrees with the expert parts within 2 mm. Nothing is modified.

Usage: python paper/margin_verification.py
Outputs: paper/outputs/margin_verification/{margin_check.csv, margin_zn_by_vertebra.csv,
         margin_summary.json, fig_margins_ALL_PLL.{pdf,png}}
"""
from __future__ import annotations

import json
import re
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
from data_io import RAW_ROOT, SUBJECT_CODE, landmark_dir, read_markup, read_vtk_polydata  # noqa: E402
from metrics import table_from_csv  # noqa: E402

OUT = REPO / "paper" / "outputs" / "margin_verification"
P2 = REPO / "part2_patient_loso_14" / "outputs"
COHORT = {"Verse005": ["L1", "L2", "L3", "L4", "L5"], "Verse561": ["L1", "L2", "L3", "L4"],
          "Verse581": ["L1", "L2", "L3", "L4", "L5"]}
CLASSES = ["ALL", "PLL"]
WALL = {"ALL": "anterior", "PLL": "posterior"}
SUP, INF, REGION, PRED, LIGHT = "#0072B2", "#D55E00", "#9a9a9a", "#111111", "#d9d9d9"
fs.VIEWS.setdefault("anterior", ((0, -1, 0), (0, 0, 1)))
fs.VIEWS.setdefault("posterior", ((0, 1, 0), (0, 0, 1)))


def parts(patient, level, lig):
    """Expert points of a class split into superior (_S) and inferior (_I) parts."""
    d = landmark_dir(patient, level)
    out = {"S": [], "I": []}
    for p in sorted(d.iterdir()):
        m = re.match(rf".*_{lig}_([SI])(_\d+)?\.mrk\.json$", p.name)
        if m:
            pts = read_markup(p)
            if pts is not None:
                out[m.group(1)].append(pts)
    return {k: np.vstack(v) if v else np.zeros((0, 3)) for k, v in out.items()}


def thresholds(held):
    train = {(s, l): {"verts": cases.mesh(s, l, "vtk"), "gt": cases.landmarks(s, l)}
             for s in COHORT if s != held for l in COHORT[s]}
    return gcz.build_thresholds(train)


def split_bands(held):
    """Counterfactual: mu +/- 2 sigma of z_n for the superior and inferior parts separately."""
    out = {}
    for lig in CLASSES:
        zs = {"S": [], "I": []}
        for s in COHORT:
            if s == held:
                continue
            for lv in COHORT[s]:
                v = np.asarray(cases.mesh(s, lv, "vtk"), float)
                for k, pts in parts(s, lv, lig).items():
                    if len(pts):
                        zs[k].append(gcz.aabb_norm(pts, v)[:, 2])
        out[lig] = {}
        for k in "SI":
            z = np.concatenate(zs[k])
            mu, sd = z.mean(), z.std(ddof=1)
            out[lig][k] = (max(0.0, mu - gcz.K_SIGMA * sd), min(1.0, mu + gcz.K_SIGMA * sd))
    return out


def region_mask(v, lig, thr):
    kd = cKDTree(v)
    d, idx = kd.query(np.asarray(gcz.candidate_vertices(v, lig, thr, third_fallback=False), float))
    m = np.zeros(len(v), bool)
    m[idx[d < 1e-6]] = True
    return m


def vertex_area(v, f):
    tri = v[f]
    a = 0.5 * np.linalg.norm(np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0]), axis=1)
    va = np.zeros(len(v))
    for c in range(3):
        np.add.at(va, f[:, c], a / 3.0)
    return va


def body_walls(v):
    """Superior and inferior margins of the anterior and posterior body walls, from the mesh only."""
    xc = 0.5 * (v[:, 0].min() + v[:, 0].max())
    mid = np.abs(v[:, 0] - xc) < 8.0
    y0 = v[:, 1].min()
    ant = mid & (v[:, 1] < y0 + 5.0)
    slab = v[np.abs(v[:, 0] - xc) < 2.0, 1]
    occ = np.zeros(int(np.ceil(slab.max() - y0)) + 2, bool)
    occ[np.floor(slab - y0).astype(int)] = True
    run, y_post = 0, None
    for i, o in enumerate(occ):
        run = 0 if o else run + 1
        if run == 4:
            y_post = y0 + i - 3
            break
    post = mid & (v[:, 1] > y_post - 3.0) & (v[:, 1] <= y_post)
    return {"xc": xc, "y_front": y0, "y_post": y_post,
            "anterior": (v[ant, 2].min(), v[ant, 2].max()), "posterior": (v[post, 2].min(), v[post, 2].max())}


def middle(z, lo_hi):
    s = (np.asarray(z, float) - lo_hi[0]) / (lo_hi[1] - lo_hi[0])
    return (s > 0.25) & (s < 0.75)


def check():
    preds = {m: table_from_csv(P2 / f"predictions_{m}.csv") for m in ("gcz", "gczrl")}
    zrows, rows = [], []
    for s in COHORT:
        for lv in COHORT[s]:
            v = np.asarray(cases.mesh(s, lv, "vtk"), float)
            for lig in CLASSES:
                for k, pts in parts(s, lv, lig).items():
                    zn = gcz.aabb_norm(pts, v)[:, 2]
                    assert np.allclose(zn, (pts[:, 2] - v[:, 2].min()) / np.ptp(v[:, 2]))
                    zrows.append({"subject": SUBJECT_CODE[s], "level": lv, "ligament": lig, "part": k,
                                  "n": len(pts), "z_n_mean": zn.mean(), "z_n_min": zn.min(), "z_n_max": zn.max()})
    for held in COHORT:
        thr, sb = thresholds(held), split_bands(held)
        for lv in COHORT[held]:
            v, f = read_vtk_polydata(RAW_ROOT / held / f"{held}_{lv}.vtk")
            v = np.asarray(v, float)
            va, walls, zn_v = vertex_area(v, f), body_walls(v), gcz.aabb_norm(v, v)[:, 2]
            for lig in CLASSES:
                pp = parts(held, lv, lig)
                reg = region_mask(v, lig, thr)
                m = (pp["I"][:, 2].mean(), pp["S"][:, 2].mean())
                band = sb[lig]
                in_split = lambda zn: ((zn >= band["S"][0]) & (zn <= band["S"][1])) | (  # noqa: E731
                    (zn >= band["I"][0]) & (zn <= band["I"][1]))
                row = {"case": f"{SUBJECT_CODE[held]} {lv}", "ligament": lig,
                       "n_sup": len(pp["S"]), "n_inf": len(pp["I"]),
                       "sup_above_inf": bool(pp["S"][:, 2].mean() > pp["I"][:, 2].mean()),
                       "sup_minus_inf_mm": m[1] - m[0],
                       "mesh_wall_height_mm": walls[WALL[lig]][1] - walls[WALL[lig]][0],
                       "band_z_n": f"{thr[lig]['z_n']['lo']:.2f}-{thr[lig]['z_n']['hi']:.2f}",
                       "split_band_S": f"{band['S'][0]:.2f}-{band['S'][1]:.2f}",
                       "split_band_I": f"{band['I'][0]:.2f}-{band['I'][1]:.2f}",
                       "region_area_middle_half_pct": 100 * va[reg & middle(v[:, 2], m)].sum() / va[reg].sum(),
                       "region_area_kept_by_split_bands_pct": 100 * va[reg & in_split(zn_v)].sum() / va[reg].sum(),
                       "expert_middle_half_pct": 100 * middle(np.vstack([pp["S"], pp["I"]])[:, 2], m).mean()}
                for name, tab in preds.items():
                    p = np.atleast_2d(tab[(held, lv)][lig])
                    row[f"{name}_points"] = len(p)
                    row[f"{name}_middle_half"] = int(middle(p[:, 2], m).sum())
                    row[f"{name}_outside_split_bands"] = int((~in_split(gcz.aabb_norm(p, v)[:, 2])).sum())
                rows.append(row)
    return pd.DataFrame(rows), pd.DataFrame(zrows)


def section_panel(ax, v, walls, regs, pp, rl):
    """Mid-sagittal section (y, z in mm): mesh slab, region vertices, body margins, and the points
    within 6 mm of the mid-sagittal plane."""
    slab = np.abs(v[:, 0] - walls["xc"]) < 1.5
    near = lambda q: np.atleast_2d(q)[np.abs(np.atleast_2d(q)[:, 0] - walls["xc"]) < 6.0]  # noqa: E731
    ax.scatter(v[slab, 1], v[slab, 2], s=1.2, color=LIGHT, lw=0, rasterized=True)
    for lig in CLASSES:
        r = slab & regs[lig]
        ax.scatter(v[r, 1], v[r, 2], s=2.2, color=REGION, lw=0, rasterized=True)
    for lig in CLASSES:
        lo, hi = walls[WALL[lig]]
        y = walls["y_front"] if lig == "ALL" else walls["y_post"]
        for z in (lo, hi):
            ax.plot([y - 7, y + 7], [z, z], color=fs.INK, lw=0.7, ls=(0, (3, 2)))
        for k, col in (("S", SUP), ("I", INF)):
            q = near(pp[lig][k])
            ax.scatter(q[:, 1], q[:, 2], s=14, color=col, lw=0, zorder=3)
        q = near(rl[lig])
        ax.scatter(q[:, 1], q[:, 2], s=7, color=PRED, lw=0, zorder=4)
        ax.text(y, hi + 3.2, lig, ha="center", fontsize=7, fontweight="bold")
    ax.set_aspect("equal")
    ax.set_xlabel("Antero-posterior coordinate y (mm), anterior to the left")
    ax.set_ylabel("Superior-inferior coordinate z (mm)")
    ax.set_xlim(walls["y_front"] - 10, walls["y_post"] + 9)
    zl = walls["anterior"][0] - 8, walls["anterior"][1] + 8
    ax.set_ylim(*zl)


def zn_panel(ax, zdf, lig, thr_by_fold, held):
    sub = zdf[zdf.ligament == lig].reset_index(drop=True)
    keys = list(dict.fromkeys(zip(sub.subject, sub.level)))
    for i, (sj, lv) in enumerate(keys):
        for k, col in (("S", SUP), ("I", INF)):
            r = sub[(sub.subject == sj) & (sub.level == lv) & (sub.part == k)].iloc[0]
            ax.plot([r.z_n_min, r.z_n_max], [i, i], color=col, lw=1.6, solid_capstyle="butt")
            ax.plot(r.z_n_mean, i, "o", color=col, ms=3.2, mec="white", mew=0.4)
    b = thr_by_fold[held][lig]["z_n"]
    ax.axvspan(b["lo"], b["hi"], color=REGION, alpha=0.16, lw=0, zorder=0)
    ax.set_yticks(range(len(keys)), [f"{sj} {lv}" for sj, lv in keys], fontsize=5.8)
    ax.set_ylim(len(keys) - 0.4, -0.6)
    held_code = SUBJECT_CODE[held]
    for t in ax.get_yticklabels():
        if t.get_text().startswith(held_code):
            t.set_fontweight("bold")
    ax.set_xlim(-0.02, 1.02)
    ax.set_xlabel(r"Normalized height $z_n$ (0 lowest, 1 highest vertex)")
    ax.text(0.5, -0.75, rf"pooled band $\mu_L\pm2\sigma_L$ (fold S03): {b['lo']:.2f}–{b['hi']:.2f}", ha="center",
            va="bottom", fontsize=6.2, color=fs.MUTED)
    ax.tick_params(axis="y", length=0)
    for sp in ("left", "right", "top"):
        ax.spines[sp].set_visible(False)


def body_only(v, f, y_cut):
    keep = np.all(v[f][:, :, 1] < y_cut, axis=1)
    return f[keep]


def figure(zdf):
    fs.setup()
    plt = fs.plt
    key, held = ("Verse581", "L4"), "Verse581"
    v, f = read_vtk_polydata(RAW_ROOT / key[0] / f"{key[0]}_{key[1]}.vtk")
    v, f = np.asarray(v, float), np.asarray(f)
    thr_by_fold = {h: thresholds(h) for h in COHORT}
    thr = thr_by_fold[held]
    rl = table_from_csv(P2 / "predictions_gczrl.csv")[key]
    walls = body_walls(v)
    regs = {lig: region_mask(v, lig, thr) for lig in CLASSES}
    pp = {lig: parts(*key, lig) for lig in CLASSES}
    fig = plt.figure(figsize=(7.2, 5.0))
    outer = fig.add_gridspec(2, 1, height_ratios=[1, 1.05], hspace=0.36)
    gs = outer[0].subgridspec(1, 3, width_ratios=[1.25, 1, 1], wspace=0.18)
    gb = outer[1].subgridspec(1, 2, wspace=0.22)

    ax = fig.add_subplot(gs[0, 0])
    section_panel(ax, v, walls, regs, pp, rl)
    ax.set_title("(a) Mid-sagittal section, S03 L4", loc="left", fontsize=7.5, fontweight="bold", pad=3)

    def margin_lines(lig, x_from, x_to, y):
        lo, hi = walls[WALL[lig]]
        return {"kind": "lines", "color": fs.INK, "radius": 0.22, "pairs": [((x_from, y, lo), (x_to, y, lo)),
                                                                          ((x_from, y, hi), (x_to, y, hi))]}

    rgb = fs.region_colors(v, {"ALL": regs["ALL"]}, order=("ALL",), mix=0.0)
    rgb[regs["ALL"]] = (np.array(fs.matplotlib.colors.to_rgb(REGION)) * 255).astype(np.uint8)
    lay = [{"kind": "spheres", "points": pp["ALL"]["S"], "color": SUP, "size": 20},
           {"kind": "spheres", "points": pp["ALL"]["I"], "color": INF, "size": 20},
           {"kind": "spheres", "points": rl["ALL"], "color": PRED, "size": 13},
           margin_lines("ALL", walls["xc"] - 26, walls["xc"] + 26, walls["y_front"] - 3)]
    ax = fig.add_subplot(gs[0, 1])
    fs.show(ax, fs.render(v, f, "anterior", lay, zoom=1.3, vertex_rgb=rgb))
    ax.set_title("(b) ALL, anterior view", loc="left", fontsize=7.5, fontweight="bold", pad=3)

    fb = body_only(v, f, walls["y_post"] + 0.5)
    rgb = fs.region_colors(v, {"PLL": regs["PLL"]}, order=("PLL",), mix=0.0)
    rgb[regs["PLL"]] = (np.array(fs.matplotlib.colors.to_rgb(REGION)) * 255).astype(np.uint8)
    lay = [{"kind": "spheres", "points": pp["PLL"]["S"], "color": SUP, "size": 20},
           {"kind": "spheres", "points": pp["PLL"]["I"], "color": INF, "size": 20},
           {"kind": "spheres", "points": rl["PLL"], "color": PRED, "size": 13},
           margin_lines("PLL", walls["xc"] - 20, walls["xc"] + 20, walls["y_post"] + 3)]
    ax = fig.add_subplot(gs[0, 2])
    fs.show(ax, fs.render(v, fb, "posterior", lay, zoom=1.3, vertex_rgb=rgb))
    ax.set_title("(c) PLL, posterior view of the body", loc="left", fontsize=7.5, fontweight="bold", pad=3)
    ax.text(0.5, -0.04, "posterior arch removed for display", transform=ax.transAxes, ha="center", va="top",
            fontsize=6.2, color=fs.MUTED)

    for c, lig in enumerate(CLASSES):
        ax = fig.add_subplot(gb[0, c])
        zn_panel(ax, zdf, lig, thr_by_fold, held)
        ax.set_title(f"({'de'[c]}) {lig}: height of the expert parts in each vertebra", loc="left", fontsize=7.5,
                     fontweight="bold", pad=12)
    items = [("Expert, superior part", SUP, "o", 6), ("Expert, inferior part", INF, "o", 6),
             ("GCZ-RL landmark", PRED, "o", 4), ("Candidate region (mask)", REGION, "s", 7)]
    handles = fs.legend_handles(items) + [fs.matplotlib.lines.Line2D([], [], color=fs.INK, lw=0.8, ls=(0, (3, 2)),
                                                                    label="Body margin (mesh)")]
    fig.legend(handles=handles, loc="upper center", ncol=5, bbox_to_anchor=(0.5, 1.0), frameon=False)
    fig.subplots_adjust(left=0.08, right=0.99, top=0.91, bottom=0.08)
    return fs.save(fig, OUT, "fig_margins_ALL_PLL", vector=False)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    df, zdf = check()
    df.round(2).to_csv(OUT / "margin_check.csv", index=False)
    zdf.to_csv(OUT / "margin_zn_by_vertebra.csv", index=False)          # full precision (checked against the text)
    summ = {"normalization_matches_direct_computation": True}
    for lig in CLASSES:
        x = df[df.ligament == lig]
        summ[lig] = {"vertebrae": len(x), "superior_above_inferior_all": bool(x.sup_above_inf.all()),
                     "sup_minus_inf_mm_median": round(float(x.sup_minus_inf_mm.median()), 1),
                     "sup_minus_inf_mm_min": round(float(x.sup_minus_inf_mm.min()), 1),
                     "region_area_middle_half_pct_median": round(float(x.region_area_middle_half_pct.median()), 1),
                     "region_area_kept_by_split_bands_pct_median": round(float(x.region_area_kept_by_split_bands_pct.median()), 1),
                     "expert_middle_half_pct_mean": round(float(x.expert_middle_half_pct.mean()), 1)}
        for name in ("gcz", "gczrl"):
            summ[lig][f"{name}_middle_half"] = f"{int(x[f'{name}_middle_half'].sum())}/{int(x[f'{name}_points'].sum())}"
            summ[lig][f"{name}_outside_split_bands"] = f"{int(x[f'{name}_outside_split_bands'].sum())}/{int(x[f'{name}_points'].sum())}"
    v = np.asarray(cases.mesh("Verse581", "L4", "vtk"), float)
    walls = body_walls(v)
    summ["figure_case_mesh_margin_minus_expert_mm"] = {
        lig: {"superior": round(float(walls[WALL[lig]][1] - parts("Verse581", "L4", lig)["S"][:, 2].mean()), 1),
              "inferior": round(float(walls[WALL[lig]][0] - parts("Verse581", "L4", lig)["I"][:, 2].mean()), 1)}
        for lig in CLASSES}
    files = figure(zdf)
    summ["figure"] = [p.name for p in files]
    (OUT / "margin_summary.json").write_text(json.dumps(summ, indent=2))
    print(json.dumps(summ, indent=1))
    print(df[["case", "ligament", "band_z_n", "split_band_S", "split_band_I", "sup_minus_inf_mm",
              "mesh_wall_height_mm", "region_area_middle_half_pct", "region_area_kept_by_split_bands_pct",
              "gcz_middle_half", "gczrl_middle_half", "gczrl_outside_split_bands"]].round(2).to_string(index=False))


if __name__ == "__main__":
    main()
