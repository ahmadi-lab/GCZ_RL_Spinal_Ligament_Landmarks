"""Analyses added in the revision (methodological comments of the co-authors).

All analyses use the primary patient-level LOSO evaluation (14 lumbar vertebrae, three folds).
They read the stored ICP geometry (part2_patient_loso_14/checkpoints) and the fold models that
part2_patient_loso_14/run_part2.py writes to outputs/models; nothing is tuned on held-out data.
Held-out annotations are read only to score the outputs, as in the main evaluation.

  A  component ablation: GCZ, + raw residual, + bounded residual, + projection (alpha = 1),
     complete GCZ-RL (class-specific alpha, second projection, side constraint); the complete
     configuration is checked against the saved GCZ-RL predictions.  In addition, the whole
     nested procedure is re-run without the post hoc side constraint.
  B  candidate-region coverage and compactness: containment of the reference landmarks,
     region vertex count, vertex fraction and surface-area fraction (triangle areas).
  C  residual interpretability: agreement of predicted and true residual vectors, grouped
     permutation importance and impurity importance of the fold regressors.
  D  atlas selection and ICP: selected atlas, level, residual, per-vertebra errors, and the
     atlas changes of the second computing environment (paper/outputs/sensitivity_fresh_icp).

Outputs: paper/outputs/revision/*.csv and revision_numbers.json.
Usage:   python paper/revision_analyses.py [part ...]   (after python part2_patient_loso_14/run_part2.py)
         parts: atlas compactness ablation interpretability no_side_constraint (default: all)
"""
from __future__ import annotations

import json
import pickle
import re
import sys
import tempfile
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
import pyvista as pv
from scipy.spatial import cKDTree
from scipy.stats import spearmanr, wilcoxon

REPO = Path(__file__).resolve().parents[1]
P2 = REPO / "part2_patient_loso_14"
OUT = REPO / "paper" / "outputs" / "revision"
sys.path.insert(0, str(REPO / "code" / "ligament"))
sys.path.insert(0, str(P2))
import cases  # noqa: E402
import gcz  # noqa: E402
import gczrl  # noqa: E402
import run_part2  # noqa: E402  (sets gczrl.SIDE_LOCKED = {"ITL", "CL"} as in the primary run)
from data_io import LIGAMENTS, RAW_ROOT, SUBJECT_CODE, read_vtk_polydata  # noqa: E402
from metrics import evaluate, table_from_csv  # noqa: E402

COHORT = run_part2.COHORT
SIDE_LOCKED = set(gczrl.SIDE_LOCKED)
KEYS = [(s, lv) for s in COHORT for lv in COHORT[s]]
FEATURES = (["x_n", "y_n", "z_n", "|x_n-0.5|", "y_n (repeat)", "z_n (repeat)",
             "distance to region centroid", "distance to region box boundary", "position in y_n band",
             "region centroid x_n", "region centroid y_n", "region centroid z_n",
             "region spread x", "region spread y", "region spread z", "region size",
             "region shape 1", "region shape 2", "region shape 3", "vertebral level"]
            + [f"class {lig}" for lig in LIGAMENTS]
            + ["surface normal x", "surface normal y", "surface normal z", "surface variation",
               "signed lateral position", "distance to class prior",
               "direction to prior x", "direction to prior y", "direction to prior z"])
GROUPS = {"Normalized position": [0, 1, 2, 3, 4, 5, 31],
          "Candidate-region relation": list(range(6, 19)),
          "Anatomical context (level, class)": list(range(19, 27)),
          "Local surface geometry": [27, 28, 29, 30],
          "Offset from class prior": [32, 33, 34, 35]}
STAGES = ["GCZ", "GCZ + projection only (no residual)", "GCZ + raw residual", "GCZ + bounded residual",
          "GCZ + bounded residual + projection", "GCZ-RL (complete)"]
N_PERM = 30


def r(x, d=3):
    return None if x is None or not np.isfinite(x) else round(float(x), d)


def fold_tag(subj):
    return f"fold_held_{SUBJECT_CODE[subj]}"


def fold_state(subj):
    geo = pickle.loads((P2 / "checkpoints" / f"{fold_tag(subj)}_gcz.pkl").read_bytes())
    m = joblib.load(P2 / "outputs" / "models" / f"{fold_tag(subj)}_extra_trees.joblib")
    ctx = gczrl.Context(thresholds=m["thresholds"], bench_norm=m["benchmark_prior"])
    return geo, m, ctx


def mesh_faces(key):
    v, f = read_vtk_polydata(RAW_ROOT / key[0] / f"{key[0]}_{key[1]}.vtk")
    assert np.array_equal(v, cases.mesh(*key, "vtk"))
    return v, f


def gt_all():
    return {k: cases.landmarks(*k) for k in KEYS}


def score(pred, label, gt):
    df = evaluate(pred, gt, label, subject_of=lambda k: SUBJECT_CODE[k[0]])
    ok = df[~df.missing_prediction]
    g = ok.groupby("ligament")
    return df, {"macro_mae_mm": r(g.mae_mm.mean().reindex(LIGAMENTS).mean()),
                "pass5_case_avg_pct": r(g.pass5.mean().reindex(LIGAMENTS).mean(), 1),
                "ref_to_pred_mm": r(ok.ref_to_pred_mean_mm.mean()),
                "pred_to_ref_mm": r(ok.pred_to_ref_mean_mm.mean()),
                "ref_within_5mm_pct": r(df.ref_pass5_all.mean(), 1),
                **{f"MAE_{lig}": r(g.mae_mm.mean()[lig]) for lig in LIGAMENTS}}


def macro_by_case(df):
    return df[~df.missing_prediction].set_index(["case", "ligament"]).mae_mm


# ------------------------------------------------------------------------ A. ablation
def ablation(res, gt):
    stages = {s: {} for s in STAGES}
    admissible = {s: [0, 0] for s in STAGES}
    surf_dist = {s: [] for s in STAGES}
    clipped = [0, 0]
    for subj in COHORT:
        geo, m, ctx = fold_state(subj)
        for key in [k for k in KEYS if k[0] == subj]:
            v, f = mesh_faces(key)
            surf = pv.PolyData(v, np.hstack([np.full((len(f), 1), 3), f]).ravel())
            base = geo["test_gcz"][key] or {}
            for s in STAGES:
                stages[s][key] = {}
            for lig in LIGAMENTS:
                gp = base.get(lig)
                if gp is None:
                    continue
                gp = np.atleast_2d(gp).astype(float)
                name, alpha, _ = m["selection"][lig]
                tau = float(m["correction_radius_mm"][lig])
                feats = np.array([gczrl.point_features(ctx, v, p, lig, key[1], key=key) for p in gp])
                dp = np.atleast_2d(m["models"][name].predict(feats))
                out = {"GCZ": gp,
                       "GCZ + projection only (no residual)": gczrl.project_to_zone(ctx, v, gp, lig, ref=gp),
                       "GCZ + raw residual": gp + dp,
                       "GCZ + bounded residual": gp + np.array([gczrl.clip(d, tau) for d in dp])}
                out["GCZ + bounded residual + projection"] = gczrl.project_to_zone(
                    ctx, v, out["GCZ + bounded residual"], lig, ref=gp)
                out["GCZ-RL (complete)"] = gczrl.blend(ctx, gp, out["GCZ + bounded residual + projection"],
                                                       alpha, v, lig)
                zone = {tuple(x) for x in ctx.zone(v, lig)}
                clipped[0] += int(sum(np.linalg.norm(d) > tau for d in dp))
                clipped[1] += len(dp)
                for s, pts in out.items():
                    stages[s][key][lig] = pts
                    admissible[s][0] += sum(tuple(p) in zone for p in pts)
                    admissible[s][1] += len(pts)
                    _, closest = surf.find_closest_cell(np.atleast_2d(pts), return_closest_point=True)
                    surf_dist[s] += list(np.linalg.norm(np.atleast_2d(closest) - pts, axis=1))
    saved = table_from_csv(P2 / "outputs" / "predictions_gczrl.csv")
    dev = max(float(np.abs(np.atleast_2d(stages["GCZ-RL (complete)"][k][lig]) - saved[k][lig]).max())
              for k in KEYS for lig in saved[k])
    res["ablation_complete_vs_saved_max_dev_mm"] = dev
    res["ablation_clipped_points"] = {"clipped": clipped[0], "points": clipped[1]}
    rows, dfs = [], {}
    for s in STAGES:
        df, sc = score(stages[s], s, gt)
        dfs[s] = df
        rows.append({"configuration": s, **sc,
                     "admissible_vertex_pct": r(100 * admissible[s][0] / admissible[s][1], 1),
                     "mean_distance_to_surface_mm": r(np.mean(surf_dist[s]))})
    gcz_case = macro_by_case(dfs["GCZ"])
    for row in rows:
        if row["configuration"] != "GCZ":
            c = macro_by_case(dfs[row["configuration"]])
            row["wilcoxon_p_vs_GCZ"] = float(wilcoxon(c, gcz_case).pvalue)
    tab = pd.DataFrame(rows)
    tab.to_csv(OUT / "table_ablation_components.csv", index=False)
    pd.concat(dfs.values()).to_csv(OUT / "ablation_per_case.csv", index=False)
    res["ablation"] = rows
    # complete GCZ-RL against the strongest geometry-only configuration (GCZ projected, no residual)
    res["GCZRL_vs_projection_only"] = paired(macro_by_case(dfs["GCZ-RL (complete)"]),
                                             macro_by_case(dfs["GCZ + projection only (no residual)"]))
    return stages


def paired(a, b, n_boot=5000):
    """Difference a - b of the macro-averaged MAE with a vertebra-level cluster bootstrap and the
    Wilcoxon test over vertebra-class cases (the conventions of paper/make_paper_numbers.py)."""
    pv = pd.concat([a.rename("a"), b.rename("b")], axis=1).dropna()

    def macro(s):
        return s.groupby(level="ligament").mean().mean()
    rng = np.random.default_rng(0)
    keys = pv.index.get_level_values("case").unique().to_numpy()
    blocks = {k: pv.xs(k, level="case", drop_level=False) for k in keys}
    boot = []
    for _ in range(n_boot):
        s = pd.concat([blocks[k] for k in rng.choice(keys, len(keys), replace=True)])
        boot.append(macro(s.a) - macro(s.b))
    d = pv.a - pv.b
    return {"difference_mm": r(macro(pv.a) - macro(pv.b)),
            "bootstrap_ci95_mm": [r(np.percentile(boot, 2.5)), r(np.percentile(boot, 97.5))],
            "wilcoxon_p": float(wilcoxon(pv.a, pv.b).pvalue), "cases": int(len(pv)),
            "lower": int((d < 0).sum()), "higher": int((d > 0).sum()), "ties": int((d == 0).sum())}


def no_side_constraint(res, gt):
    """Whole nested procedure (inner selection included) re-run without the side constraint."""
    with tempfile.TemporaryDirectory() as tmp:
        old_out, run_part2.OUT = run_part2.OUT, Path(tmp)
        (Path(tmp) / "models").mkdir()
        gczrl.SIDE_LOCKED = set()
        try:
            pred, sel = {}, {}
            for subj in COHORT:
                _, p, info, _ = run_part2.run_fold(subj, fresh=False)
                pred.update(p)
                sel[info["held_out_subject"]] = {l: [v["extra_trees"], v["alpha"]] for l, v in info["selection"].items()}
        finally:
            gczrl.SIDE_LOCKED = SIDE_LOCKED
            run_part2.OUT = old_out
    _, sc = score(pred, "GCZ-RL without side constraint", gt)
    same_side = 0
    for k in KEYS:
        p = (pred.get(k) or {}).get("ITL")
        if p is not None and len(p) == 2:
            cx = float(cases.mesh(*k, "vtk")[:, 0].mean())
            same_side += int((p[0, 0] >= cx) == (p[1, 0] >= cx))
    res["no_side_constraint"] = {**sc, "ITL_cases_both_points_one_side": same_side, "selection": sel}
    pd.DataFrame([{"configuration": "GCZ-RL without side constraint (nested re-run)", **sc,
                   "ITL_cases_both_points_one_side": same_side}]).to_csv(
        OUT / "table_ablation_no_side_constraint.csv", index=False)


# ------------------------------------------------------------- B. region compactness
def compactness(res, gt):
    rows = []
    for subj in COHORT:
        _, m, ctx = fold_state(subj)
        thr = m["thresholds"]
        for key in [k for k in KEYS if k[0] == subj]:
            v, f = mesh_faces(key)
            tri = v[f]
            area = 0.5 * np.linalg.norm(np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0]), axis=1)
            vert_area = np.zeros(len(v))
            for c in range(3):
                np.add.at(vert_area, f[:, c], area / 3.0)
            kd = cKDTree(v)
            for lig in LIGAMENTS:
                g = (gt[key] or {}).get(lig)
                if g is None or len(g) == 0:
                    continue
                for zone_kind, third in (("selection", True), ("projection", False)):
                    zset = {tuple(x) for x in gcz.candidate_vertices(v, lig, thr, third_fallback=third)}
                    mask = np.fromiter((tuple(x) in zset for x in v), bool, len(v))
                    inside = mask[kd.query(np.atleast_2d(g))[1]]
                    rows.append({"held_out": SUBJECT_CODE[subj], "case": f"{SUBJECT_CODE[subj]}_{key[1]}",
                                 "ligament": lig, "zone": zone_kind, "n_ref": len(g), "n_inside": int(inside.sum()),
                                 "mesh_vertices": len(v), "zone_vertices": int(mask.sum()),
                                 "vertex_fraction_pct": 100 * mask.mean(),
                                 "mesh_area_mm2": float(area.sum()),
                                 "area_fraction_pct": 100 * vert_area[mask].sum() / vert_area.sum()})
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "region_compactness_per_vertebra.csv", index=False)

    def q(x):
        return f"{np.median(x):.1f} ({np.percentile(x, 25):.1f}-{np.percentile(x, 75):.1f})"
    tab = []
    for lig in LIGAMENTS:
        s, p = df[(df.ligament == lig) & (df.zone == "selection")], df[(df.ligament == lig) & (df.zone == "projection")]
        tab.append({"ligament": lig, "vertebrae": len(p), "n_ref": int(p.n_ref.sum()),
                    "inside_selection_pct": r(100 * s.n_inside.sum() / s.n_ref.sum(), 1),
                    "inside_projection_pct": r(100 * p.n_inside.sum() / p.n_ref.sum(), 1),
                    "projection_vertices_median": float(np.median(p.zone_vertices)),
                    "projection_vertex_fraction_pct": q(p.vertex_fraction_pct),
                    "projection_area_fraction_pct": q(p.area_fraction_pct),
                    "projection_area_fraction_median": r(np.median(p.area_fraction_pct), 1),
                    "projection_area_fraction_by_subject": {h: r(np.median(x.area_fraction_pct), 1)
                                                            for h, x in p.groupby("held_out")},
                    "selection_area_fraction_median": r(np.median(s.area_fraction_pct), 1)})
    pd.DataFrame(tab).to_csv(OUT / "table_region_compactness.csv", index=False)
    s = df[df.zone == "selection"]
    p = df[df.zone == "projection"]
    res["compactness"] = {"by_class": tab,
                          "inside_selection_pct_all": r(100 * s.n_inside.sum() / s.n_ref.sum(), 1),
                          "inside_projection_pct_all": r(100 * p.n_inside.sum() / p.n_ref.sum(), 1),
                          "mesh_vertices_median": float(np.median(p.drop_duplicates("case").mesh_vertices)),
                          "mesh_area_median_mm2": r(np.median(p.drop_duplicates("case").mesh_area_mm2), 0),
                          "area_fraction_median_all_cases": r(np.median(p.area_fraction_pct), 1)}


# ------------------------------------------------------- C. residual interpretability
def interpretability(res, gt):
    rows, imp_perm, imp_mdi = [], [], []
    rng = np.random.default_rng(0)
    for subj in COHORT:
        geo, m, ctx = fold_state(subj)
        X, Y, names = [], [], []
        for key in [k for k in KEYS if k[0] == subj]:
            v = cases.mesh(*key, "vtk")
            for lig in LIGAMENTS:
                gp, g = (geo["test_gcz"][key] or {}).get(lig), (gt[key] or {}).get(lig)
                if gp is None or g is None or len(g) == 0:
                    continue
                gp = np.atleast_2d(gp).astype(float)
                true = gczrl.residuals(gp, g)
                name, alpha, _ = m["selection"][lig]
                tau = float(m["correction_radius_mm"][lig])
                feats = np.array([gczrl.point_features(ctx, v, p, lig, key[1], key=key) for p in gp])
                dp = np.atleast_2d(m["models"][name].predict(feats))
                for i in range(len(gp)):
                    c = gczrl.clip(dp[i], tau)
                    nt, npd = np.linalg.norm(true[i]), np.linalg.norm(dp[i])
                    rows.append({"held_out": SUBJECT_CODE[subj], "case": f"{SUBJECT_CODE[subj]}_{key[1]}",
                                 "ligament": lig, "point": i, "alpha": alpha, "model": name,
                                 "true_mm": nt, "pred_mm": npd, "pred_clipped_mm": float(np.linalg.norm(c)),
                                 "cosine": float(true[i] @ dp[i] / (nt * npd)) if nt > 1e-9 and npd > 1e-9 else np.nan,
                                 "dist_before_mm": nt, "dist_after_clipped_mm": float(np.linalg.norm(true[i] - c)),
                                 **{f"true_{a}": true[i][j] for j, a in enumerate("xyz")},
                                 **{f"pred_{a}": dp[i][j] for j, a in enumerate("xyz")},
                                 **{f"gcz_{a}": gp[i][j] for j, a in enumerate("xyz")}})
                    X.append(feats[i])
                    Y.append(true[i])
                    names.append(name)
        X, Y, names = np.array(X), np.array(Y), np.array(names)

        def err(Xm):
            pred = np.zeros_like(Y)
            for nm in np.unique(names):
                idx = names == nm
                pred[idx] = m["models"][nm].predict(Xm[idx])
            return float(np.linalg.norm(pred - Y, axis=1).mean())
        base = err(X)
        for grp, cols in GROUPS.items():
            inc = []
            for _ in range(N_PERM):
                Xp = X.copy()
                Xp[:, cols] = X[rng.permutation(len(X))][:, cols]
                inc.append(err(Xp) - base)
            imp_perm.append({"held_out": SUBJECT_CODE[subj], "group": grp, "baseline_error_mm": base,
                             "increase_mm": float(np.mean(inc))})
        for nm, model in m["models"].items():
            fi = model.feature_importances_
            imp_mdi.append({"held_out": SUBJECT_CODE[subj], "model": nm,
                            **{grp: float(fi[cols].sum()) for grp, cols in GROUPS.items()},
                            **{f"f{i + 1}": float(x) for i, x in enumerate(fi)}})
    pts = pd.DataFrame(rows)
    pts.to_csv(OUT / "residual_vectors_heldout.csv", index=False)
    perm = pd.DataFrame(imp_perm)
    perm.to_csv(OUT / "feature_importance_permutation.csv", index=False)
    mdi = pd.DataFrame(imp_mdi)
    mdi.to_csv(OUT / "feature_importance_impurity.csv", index=False)
    by_cls = []
    for lig in LIGAMENTS:
        x = pts[pts.ligament == lig]
        by_cls.append({"ligament": lig, "points": len(x), "true_residual_mean_mm": r(x.true_mm.mean(), 2),
                       "pred_clipped_mean_mm": r(x.pred_clipped_mm.mean(), 2),
                       "median_cosine": r(np.nanmedian(x.cosine), 2),
                       "pct_cosine_positive": r(100 * (x.cosine > 0).mean(), 1),
                       "pct_closer_after_clipped": r(100 * (x.dist_after_clipped_mm < x.dist_before_mm).mean(), 1)})
    pd.DataFrame(by_cls).to_csv(OUT / "table_residual_agreement.csv", index=False)
    g = perm.groupby("group").increase_mm
    top = mdi[[f"f{i + 1}" for i in range(36)]].mean().sort_values(ascending=False).head(10)
    res["interpretability"] = {
        "points": len(pts), "by_class": by_cls,
        "pct_cosine_positive_all": r(100 * (pts.cosine > 0).mean(), 1),
        "pct_closer_after_clipped_all": r(100 * (pts.dist_after_clipped_mm < pts.dist_before_mm).mean(), 1),
        "permutation_increase_mm": {k: {"mean": r(g.mean()[k]), "min": r(g.min()[k]), "max": r(g.max()[k])}
                                    for k in GROUPS},
        "permutation_baseline_error_mm": {h: r(x.baseline_error_mm.iloc[0]) for h, x in perm.groupby("held_out")},
        "impurity_group_share": {k: r(mdi[k].mean()) for k in GROUPS},
        "impurity_top10": [{"feature": FEATURES[int(k[1:]) - 1], "importance": r(v)} for k, v in top.items()]}


# ------------------------------------------------------------- D. atlas and ICP
def atlas(res):
    def choices(path):
        out = {}
        for f in json.loads(Path(path).read_text()):
            for tgt, s in f["test_atlas_choice"].items():
                mm = re.match(r"(\S+) \(residual ([\d.]+) mm\)", s)
                out[tgt] = (mm.group(1), float(mm.group(2)))
        return out
    stored = choices(P2 / "outputs" / "fold_details.json")
    fresh_dir = REPO / "paper" / "outputs" / "sensitivity_fresh_icp"
    fresh = choices(fresh_dir / "fold_details.json")
    vert = pd.read_csv(REPO / "paper" / "outputs" / "primary_per_vertebra_macro_mae.csv").set_index("case")
    fm = pd.read_csv(fresh_dir / "per_vertebra_per_ligament_metrics.csv")
    fvert = fm[~fm.missing_prediction].groupby(["case", "method"]).mae_mm.mean().unstack()
    src = {v: k for k, v in SUBJECT_CODE.items()}
    rows = []
    for tgt, (atl, resid) in stored.items():
        case = f"{src[tgt.split('_')[0]]}_{tgt.split('_')[1]}"
        rows.append({"target": tgt, "atlas": atl, "same_level": atl.split("_")[1] == tgt.split("_")[1],
                     "icp_residual_mm": resid, "GCZ_mm": vert.loc[case, "GCZ"], "GCZRL_mm": vert.loc[case, "GCZ-RL"],
                     "SLD_mm": vert.loc[case, "SLD"], "second_env_atlas": fresh[tgt][0],
                     "second_env_residual_mm": fresh[tgt][1], "atlas_changed": fresh[tgt][0] != atl,
                     "second_env_GCZ_mm": r(fvert.loc[case, "GCZ"]), "second_env_GCZRL_mm": r(fvert.loc[case, "GCZ-RL"])})
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "table_atlas_icp.csv", index=False)
    same, cross = df[df.same_level], df[~df.same_level]
    rho_g = spearmanr(df.icp_residual_mm, df.GCZ_mm)
    rho_r = spearmanr(df.icp_residual_mm, df.GCZRL_mm)
    res["atlas"] = {"same_level": int(len(same)), "cross_level": int(len(cross)),
                    "same_level_GCZ_mm": r(same.GCZ_mm.mean()), "cross_level_GCZ_mm": r(cross.GCZ_mm.mean()),
                    "same_level_GCZRL_mm": r(same.GCZRL_mm.mean()), "cross_level_GCZRL_mm": r(cross.GCZRL_mm.mean()),
                    "spearman_residual_vs_GCZ": [r(rho_g.statistic), r(rho_g.pvalue)],
                    "spearman_residual_vs_GCZRL": [r(rho_r.statistic), r(rho_r.pvalue)],
                    "changed_in_second_env": df[df.atlas_changed].target.tolist(),
                    "GCZRL_improved_over_GCZ_vertebrae": int((df.GCZRL_mm < df.GCZ_mm).sum())}


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    parts = sys.argv[1:] or ["atlas", "compactness", "ablation", "interpretability", "no_side_constraint"]
    path = OUT / "revision_numbers.json"
    res = json.loads(path.read_text()) if path.exists() else {}
    gt = gt_all()
    for name in parts:
        fn = globals()[name]
        fn(res) if name == "atlas" else fn(res, gt)
    path.write_text(json.dumps(res, indent=2, default=float), encoding="utf-8")
    ok = res.get("ablation_complete_vs_saved_max_dev_mm", 0.0) < 1e-9
    print("complete configuration reproduces the saved GCZ-RL predictions:", ok)
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
