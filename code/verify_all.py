"""Independent checks of both parts.  Run after run_part1.py and run_part2.py.

Writes VERIFICATION_REPORT.json in the package root; exit code 1 if any check fails.
Metrics are recomputed here with plain numpy/scipy from the written prediction files.
"""
import json
import sys
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment
from sklearn.ensemble import ExtraTreesRegressor

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code" / "ligament"))
from data_io import LIGAMENTS, SOURCE_ID, load_landmarks  # noqa: E402

P1, P2 = ROOT / "part1_dissertation_57" / "outputs", ROOT / "part2_patient_loso_14" / "outputs"
RESULTS = []


def check(name, ok, detail=""):
    RESULTS.append({"check": name, "passed": bool(ok), "detail": detail})
    print(f"{'PASS' if ok else 'FAIL'}  {name}  {detail}")


def load(path):
    d = pd.read_csv(path)
    return {(c, l): g.sort_values("point_index")[["x", "y", "z"]].to_numpy() for (c, l), g in d.groupby(["case", "ligament"])}


def recompute(pred, cases):
    mae, hit, n = {l: [] for l in LIGAMENTS}, {l: 0 for l in LIGAMENTS}, {l: 0 for l in LIGAMENTS}
    for case in cases:
        code, lv = case.split("_")
        gt = load_landmarks(SOURCE_ID[code], lv)
        for lig in LIGAMENTS:
            p = pred.get((case, lig))
            if p is None or gt.get(lig) is None:
                continue
            d = np.linalg.norm(p[:, None] - gt[lig][None], axis=2)
            r, c = linear_sum_assignment(d)
            mae[lig].append(d[r, c].mean())
            hit[lig] += int((d[r, c] <= 5).sum())
            n[lig] += len(r)
    return float(np.mean([np.mean(mae[l]) for l in LIGAMENTS])), float(np.mean([hit[l] / n[l] * 100 for l in LIGAMENTS]))


def main():
    # ---------------- Part 1
    v1 = json.load(open(P1 / "verification_part1.json"))
    check("P1 Table 4.4 equals the dissertation (all MAE and Pass@5 values)", v1["table_4_4_matches_dissertation"])
    ag = v1["regenerated_vs_archived_predictions"]
    for m in ("GCZ", "GCZ-RL"):
        a = ag[m]
        check(f"P1 {m} predictions identical to the dissertation predictions",
              a["identical_sets_tol_1um"] == a["case_ligament_sets_compared"],
              f"{a['identical_sets_tol_1um']}/{a['case_ligament_sets_compared']} point sets")
    t44 = pd.read_csv(P1 / "table_4_4_heldout_results.csv")
    held = sorted({c for c, _ in load(P1 / "predictions_gczrl_extra_trees_heldout37.csv")
                   if load_landmarks(SOURCE_ID[c.split('_')[0]], c.split('_')[1]) is not None})
    mae, ps = recompute(load(P1 / "predictions_gczrl_extra_trees_heldout37.csv"), held)
    row = t44[(t44.method == "GCZ-RL") & (t44.ligament == "Mean")].iloc[0]
    check("P1 GCZ-RL metrics recomputed from coordinates", abs(mae - row.MAE_mm) < 5e-4 and abs(ps - row.Pass5_percent) < 5e-2,
          f"{len(held)} vertebrae, MAE {mae:.3f}, Pass@5 {ps:.1f}")
    final = joblib.load(P1 / "models" / "gczrl_extra_trees_final.joblib")
    check("P1 final models are Extra Trees", all(isinstance(m, ExtraTreesRegressor) for m in final["models"].values()),
          ", ".join(sorted(final["models"])))

    # ---------------- Part 2
    folds = json.load(open(P2 / "fold_details.json"))
    ok = len(folds) == 3 and all(not set(f["training_vertebrae"]) & set(f["test_vertebrae"]) and
                                 all(v.split("_")[0] != f["held_out_subject"] for v in f["training_vertebrae"])
                                 for f in folds)
    check("P2 three folds, held-out subject never in training", ok,
          "; ".join(f"{f['held_out_subject']}: train {len(f['training_vertebrae'])}, test {len(f['test_vertebrae'])}" for f in folds))
    et_ok, det = True, []
    for f in folds:
        bundle = joblib.load(P2 / "models" / f"{f['fold']}_extra_trees.joblib")
        types = {type(m).__name__ for m in bundle["models"].values()}
        et_ok &= types == {"ExtraTreesRegressor"} and all(s["extra_trees"].startswith("ET-") for s in f["selection"].values())
        det.append(f"{f['held_out_subject']}: {sorted(types)}")
    check("P2 every selected and saved model is ExtraTreesRegressor", et_ok, "; ".join(det))
    overall = pd.read_csv(P2 / "table_P2_2_overall.csv")
    cohort = sorted({c for c, _ in load(P2 / "predictions_gcz.csv")})
    worst = 0.0
    for m, f in [("GCZ", "predictions_gcz.csv"), ("GCZ-RL", "predictions_gczrl.csv"), ("SLD", "predictions_sld.csv")]:
        mae, ps = recompute(load(P2 / f), cohort)
        r = overall[overall.method == m].iloc[0]
        worst = max(worst, abs(mae - r.MAE_mm), abs(ps - r.Pass5_percent) / 100)
    check("P2 metrics recomputed from coordinates", worst < 5e-4 and len(cohort) == 14, f"14 vertebrae, max diff {worst:.1e}")
    leak = P2 / "table_P2_5_leakage_test.json"
    if leak.exists():
        lt = json.load(open(leak))
        check("P2 leakage test: corrupting held-out annotations leaves predictions unchanged", lt["predictions_unchanged"],
              f"max difference {lt['max_abs_coordinate_difference_mm']:.1e} mm over {lt['points_compared']} points")
    else:
        check("P2 leakage test: corrupting held-out annotations leaves predictions unchanged", False, "not run")
    out = {"all_passed": all(r["passed"] for r in RESULTS), "checks": RESULTS}
    (ROOT / "VERIFICATION_REPORT.json").write_text(json.dumps(out, indent=2), encoding="utf-8")
    return 0 if out["all_passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
