"""Evaluation metrics shared by every analysis in the package.

Matching rule (identical for all methods)
-----------------------------------------
For one vertebra and one ligament, let P be the predicted landmark set and G
the expert reference set.  A one-to-one Hungarian assignment minimising the
total Euclidean distance is computed between P and G.  Because |P| and |G|
generally differ, the assignment covers min(|P|, |G|) pairs.

Reported quantities
-------------------
* ``mae_mm`` / ``pass5``: mean matched distance and percentage of matched
  pairs <= 5 mm.  These are *conditional* on the matched pairs; they are the
  metrics used in the dissertation (Chapter 4).
* ``ref_coverage``: fraction of reference points that received a match
  (= n_matched / |G|).  Unmatched reference points are not penalised by the
  conditional metrics, so coverage is reported alongside them.
* ``ref_to_pred_mean_mm``: mean distance from *every* reference point to its
  nearest predicted point (penalises unmatched reference structure).
* ``pred_to_ref_mean_mm``: mean distance from every prediction to its
  nearest reference point.
* ``symmetric_nn_mm``: average of the two directed means above.
* ``ref_pass5_all``: percentage of *all* reference points whose nearest
  prediction lies within 5 mm (a coverage-aware pass rate).
Missing predictions (no points for a ligament) are counted explicitly and are
not silently dropped from coverage statistics.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment
from scipy.spatial import cKDTree

PASS_MM = 5.0


def pair_metrics(pred, gt) -> dict:
    out = {"n_pred": 0, "n_ref": 0, "n_matched": 0, "n_pass5": 0, "missing_prediction": True,
           "mae_mm": np.nan, "pass5": np.nan, "ref_coverage": 0.0,
           "ref_to_pred_mean_mm": np.nan, "pred_to_ref_mean_mm": np.nan,
           "symmetric_nn_mm": np.nan, "ref_pass5_all": np.nan}
    if gt is None or len(gt) == 0:
        return out
    g = np.atleast_2d(np.asarray(gt, float))
    out["n_ref"] = len(g)
    if pred is None or len(pred) == 0:
        out["ref_pass5_all"] = 0.0
        return out
    p = np.atleast_2d(np.asarray(pred, float))
    d = np.linalg.norm(p[:, None, :] - g[None, :, :], axis=2)
    ri, ci = linear_sum_assignment(d)
    md = d[ri, ci]
    r2p = cKDTree(p).query(g)[0]
    p2r = cKDTree(g).query(p)[0]
    out.update({
        "n_pred": len(p), "n_matched": len(md), "n_pass5": int((md <= PASS_MM).sum()),
        "missing_prediction": False,
        "mae_mm": float(md.mean()), "pass5": float((md <= PASS_MM).mean() * 100),
        "ref_coverage": float(len(md) / len(g)),
        "ref_to_pred_mean_mm": float(r2p.mean()), "pred_to_ref_mean_mm": float(p2r.mean()),
        "symmetric_nn_mm": float(0.5 * (r2p.mean() + p2r.mean())),
        "ref_pass5_all": float((r2p <= PASS_MM).mean() * 100),
    })
    return out


def evaluate(pred: dict, gt: dict, method: str, subject_of=lambda key: key[0]) -> pd.DataFrame:
    """Case-ligament table for every (case, ligament) that has reference points."""
    rows = []
    for key, gt_by_lig in gt.items():
        if gt_by_lig is None:
            continue
        for lig, g in gt_by_lig.items():
            if g is None or len(g) == 0:
                continue
            p = (pred.get(key) or {}).get(lig)
            m = pair_metrics(p, g)
            rows.append({"method": method, "subject": subject_of(key), "case": f"{key[0]}_{key[1]}",
                         "level": key[1], "ligament": lig, **m})
    return pd.DataFrame(rows)


def summarize(df: pd.DataFrame, ligaments) -> dict:
    """Dissertation convention (Chapter 4, verified against the archived tables):

    * MAE = macro average over ligaments of the mean case-level matched error;
    * Pass@5 = macro average over ligaments of the pooled matched-point pass rate
      (all matched pairs of a ligament pooled across cases).
    Case-ligament-averaged and fully pooled pass rates are reported as well.
    """
    ok = df[~df.missing_prediction]
    per_lig = ok.groupby("ligament").mae_mm.mean().reindex(ligaments)
    g = ok.groupby("ligament")
    per_lig_pass = (g.n_pass5.sum() / g.n_matched.sum() * 100).reindex(ligaments)
    return {
        "macro_mae_mm": float(per_lig.mean()),
        "pass5": float(per_lig_pass.mean()),
        "pass5_case_ligament_mean": float(ok.pass5.mean()),
        "pass5_all_points_pooled": float(ok.n_pass5.sum() / ok.n_matched.sum() * 100),
        "ref_coverage_mean": float(df.ref_coverage.mean()),
        "ref_to_pred_mean_mm": float(ok.ref_to_pred_mean_mm.mean()),
        "symmetric_nn_mm": float(ok.symmetric_nn_mm.mean()),
        "ref_pass5_all_mean": float(df.ref_pass5_all.mean()),
        "missing_case_ligaments": int(df.missing_prediction.sum()),
        "case_ligaments": int(len(df)),
        "cases": int(df.case.nunique()),
        "per_ligament_mae_mm": {k: float(v) for k, v in per_lig.items()},
        "per_ligament_pass5": {k: float(v) for k, v in per_lig_pass.items()},
    }


def table_from_csv(path, method_col_ok=True) -> dict:
    """Read a coordinate CSV in any of the archived layouts into {(pat, lv): {lig: pts}}."""
    df = pd.read_csv(path)
    if "pred_pts_json" in df:  # archived SLD layout: one JSON list per case-ligament
        import json
        out: dict = {}
        for _, r in df.iterrows():
            txt = r["pred_pts_json"]
            pts = np.asarray(json.loads(txt), float) if isinstance(txt, str) and txt.strip() else None
            out.setdefault((str(r["patient"]), str(r["level"])), {})[str(r["ligament"])] = (
                pts.reshape(-1, 3) if pts is not None and pts.size else None)
        return out
    pat = "patient_id" if "patient_id" in df else "patient"
    lv = "vertebra_id" if "vertebra_id" in df else "level"
    idx = "landmark_id" if "landmark_id" in df else "point_index"
    out: dict = {}
    for (p, l, lig), grp in df.groupby([pat, lv, "ligament"], sort=False):
        out.setdefault((str(p), str(l)), {})[str(lig)] = grp.sort_values(idx)[["x", "y", "z"]].to_numpy(float)
    return out


def to_rows(pred: dict, method: str, code=lambda p, l: f"{p}_{l}") -> pd.DataFrame:
    rows = []
    for (p, l), by_lig in pred.items():
        for lig, pts in (by_lig or {}).items():
            if pts is None:
                continue
            for i, (x, y, z) in enumerate(np.atleast_2d(pts)):
                rows.append({"patient": p, "level": l, "case": code(p, l), "ligament": lig,
                             "point_index": i, "x": float(x), "y": float(y), "z": float(z), "method": method})
    return pd.DataFrame(rows)


def coordinate_agreement(a: dict, b: dict) -> dict:
    """Compare two prediction dictionaries point-set by point-set.

    Uses Hungarian matching within each case-ligament so that point order does
    not matter.  Returns counts of identical sets and distance statistics.
    """
    keys = sorted(set((k, lig) for k, v in a.items() for lig in (v or {}) if v[lig] is not None)
                  | set((k, lig) for k, v in b.items() for lig in (v or {}) if v[lig] is not None))
    n_same = n_cmp = n_missing = 0
    dists = []
    for k, lig in keys:
        pa = (a.get(k) or {}).get(lig)
        pb = (b.get(k) or {}).get(lig)
        if pa is None or pb is None:
            n_missing += 1
            continue
        pa = np.atleast_2d(pa)
        pb = np.atleast_2d(pb)
        n_cmp += 1
        if pa.shape != pb.shape:
            dists.append(np.inf)
            continue
        d = np.linalg.norm(pa[:, None] - pb[None], axis=2)
        ri, ci = linear_sum_assignment(d)
        mx = float(d[ri, ci].max())
        dists.append(mx)
        if mx <= 1e-3:
            n_same += 1
    arr = np.array([x for x in dists if np.isfinite(x)])
    return {"case_ligament_sets_compared": n_cmp, "identical_sets_tol_1um": n_same,
            "sets_missing_in_one": n_missing,
            "sets_with_different_point_count": int(sum(1 for x in dists if not np.isfinite(x))),
            "max_point_deviation_mm": float(arr.max()) if len(arr) else float("nan"),
            "median_set_max_deviation_mm": float(np.median(arr)) if len(arr) else float("nan")}
