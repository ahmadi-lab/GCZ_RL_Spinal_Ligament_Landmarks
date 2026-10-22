"""Regenerate the numbers of the journal manuscript that are not written by run_part1/run_part2.

Reads the saved predictions, metrics and checkpoints of both parts (nothing is retrained) and
writes paper/outputs/paper_numbers.json plus CSV tables:

  * primary evaluation: Pass@5 under both conventions, paired Wilcoxon tests, vertebra-level
    cluster-bootstrap CIs, per-vertebra macro MAE, ITL errors, mesh sizes, atlas residuals
    (main text Sections 4.1 and 4.3);
  * Supplementary Tables S2-S4: point accounting, matching-independent metrics (closest point in
    both directions, symmetric, reference points within 5 mm), candidate-region feasibility;
  * development analysis: metrics of the GCZ run that GCZ-RL corrects (Section 4.4).

Every value quoted in the manuscript is compared with paper/expected_values.json; exit code 1
on any mismatch.  Usage:  python paper/make_paper_numbers.py
"""
from __future__ import annotations

import json
import pickle
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.stats import wilcoxon

REPO = Path(__file__).resolve().parents[1]
OUT = REPO / "paper" / "outputs"
sys.path.insert(0, str(REPO / "code" / "ligament"))
import cases  # noqa: E402
import gcz  # noqa: E402
from data_io import LIGAMENTS, read_manifest  # noqa: E402
from metrics import evaluate, summarize  # noqa: E402

P1 = REPO / "part1_dissertation_57"
P2 = REPO / "part2_patient_loso_14"
COHORT = {"Verse005": ["L1", "L2", "L3", "L4", "L5"], "Verse561": ["L1", "L2", "L3", "L4"],
          "Verse581": ["L1", "L2", "L3", "L4", "L5"]}
N_BOOT = 5000


def r(x, d=3):
    return round(float(x), d)


def primary(res):
    df = pd.read_csv(P2 / "outputs" / "per_vertebra_per_ligament_metrics.csv")
    ok = df[~df.missing_prediction]
    s3, t4 = [], []
    for m in ["GCZ", "GCZ-RL", "SLD"]:
        d = ok[ok.method == m]
        g = d.groupby("ligament")
        case_avg = g.pass5.mean().reindex(LIGAMENTS)
        pooled = (g.n_pass5.sum() / g.n_matched.sum() * 100).reindex(LIGAMENTS)
        mae = g.mae_mm.mean().reindex(LIGAMENTS)
        t4 += [{"method": m, "ligament": lig, "MAE_mm": r(mae[lig]), "Pass5_case_avg": r(case_avg[lig], 1),
                "Pass5_pooled": r(pooled[lig], 1)} for lig in LIGAMENTS]
        res[m] = {"macro_mae": r(g.mae_mm.mean().mean()),
                  "pass5_case_avg": r(case_avg.mean(), 1), "pass5_pooled_within_class": r(pooled.mean(), 1),
                  "n_ref": int(d.n_ref.sum()), "n_pred": int(d.n_pred.sum()),
                  "pooled_coverage_pct": r(d.n_matched.sum() / d.n_ref.sum() * 100, 1),
                  "closest_point_pred_to_ref": r(d.pred_to_ref_mean_mm.mean()),
                  "closest_point_ref_to_pred": r(d.ref_to_pred_mean_mm.mean()),
                  "symmetric_nn": r(d.symmetric_nn_mm.mean()),
                  "ref_within_5mm_pct": r(d.ref_pass5_all.mean(), 1)}
        for lig in LIGAMENTS:
            x = d[d.ligament == lig]
            s3.append({"method": m, "ligament": lig, "P_to_R": r(x.pred_to_ref_mean_mm.mean()),
                       "R_to_P": r(x.ref_to_pred_mean_mm.mean()), "symmetric": r(x.symmetric_nn_mm.mean()),
                       "ref_within_5mm_pct": r(x.ref_pass5_all.mean(), 1),
                       "n_ref": int(x.n_ref.sum()), "n_pred": int(x.n_pred.sum())})
    pd.DataFrame(s3).to_csv(OUT / "table_S2_S3_point_accounting_and_matching_independent.csv", index=False)
    pd.DataFrame(t4).to_csv(OUT / "table4_primary_per_ligament.csv", index=False)

    pv = ok.pivot_table(index=["case", "ligament"], columns="method", values="mae_mm")
    vert = pv.groupby(level="case").mean()
    res["vertebrae_improved_GCZRL_vs_GCZ"] = int((vert["GCZ-RL"] < vert["GCZ"] - 1e-9).sum())
    vert.round(3).to_csv(OUT / "primary_per_vertebra_macro_mae.csv")
    itl = pv.xs("ITL", level="ligament")["GCZ-RL"].sort_values()
    res["ITL_GCZRL_errors_mm"] = {k: r(v, 2) for k, v in itl.items()}

    def macro(sub, m):
        return sub[m].groupby(level="ligament").mean().reindex(LIGAMENTS).mean()

    rng = np.random.default_rng(0)
    keys = vert.index.to_numpy()
    blocks = {k: pv.xs(k, level="case", drop_level=False) for k in keys}
    boot = {"GCZ": [], "SLD": []}
    for _ in range(N_BOOT):
        sub = pd.concat([blocks[k] for k in rng.choice(keys, len(keys), replace=True)])
        for other in boot:
            boot[other].append(macro(sub, "GCZ-RL") - macro(sub, other))
    for other in boot:
        diff = (pv["GCZ-RL"] - pv[other]).dropna()
        res[f"GCZRL_vs_{other}"] = {
            "macro_mae_difference_mm": r(macro(pv, "GCZ-RL") - macro(pv, other)),
            "bootstrap_ci95_mm": [r(np.percentile(boot[other], 2.5)), r(np.percentile(boot[other], 97.5))],
            "wilcoxon_p": float(wilcoxon(pv["GCZ-RL"], pv[other]).pvalue),
            "cases": int(len(diff)), "GCZRL_lower": int((diff < 0).sum()), "other_lower": int((diff > 0).sum()),
            "ties": int((diff == 0).sum())}

    res["mesh_vertices_median"] = float(np.median([len(cases.mesh(p, l, "vtk")) for p, ls in COHORT.items() for l in ls]))
    folds = json.loads((P2 / "outputs" / "fold_details.json").read_text())
    resid, same_level = [], 0
    for f in folds:
        for tgt, s in f["test_atlas_choice"].items():
            resid.append(float(re.search(r"residual ([\d.]+)", s).group(1)))
            same_level += tgt.split("_")[1] == s.split()[0].split("_")[1]
    res["atlas_residual_mm"] = {"mean": r(np.mean(resid), 2), "min": r(min(resid), 2), "max": r(max(resid), 2)}
    res["atlas_same_level"] = same_level


def zones(res):
    rows = []
    for held in COHORT:
        train = {(s, l): {"verts": cases.mesh(s, l, "vtk"), "gt": cases.landmarks(s, l)}
                 for s in COHORT if s != held for l in COHORT[s]}
        thr = gcz.build_thresholds(train)
        for lv in COHORT[held]:
            v, gt = cases.mesh(held, lv, "vtk"), cases.landmarks(held, lv)
            ff = gcz.threshold_features(gcz.aabb_norm(v, v))
            nn_idx = cKDTree(v)
            for lig in LIGAMENTS:
                g = gt.get(lig)
                if g is None or len(g) == 0:
                    continue
                full = np.ones(len(v), bool)
                for f in gcz.THRESHOLD_FEATURES:
                    full &= (ff[f] >= thr[lig][f]["lo"]) & (ff[f] <= thr[lig][f]["hi"])
                sel = gcz.candidate_vertices(v, lig, thr, third_fallback=True)
                proj = gcz.candidate_vertices(v, lig, thr, third_fallback=False)
                sel_set = {tuple(x) for x in sel}
                inside = np.array([tuple(x) in sel_set for x in v[nn_idx.query(g)[1]]])
                cx = float(v[:, 0].mean())
                empty = lig in ("ITL", "CL") and ((proj[:, 0] >= cx).sum() == 0 or (proj[:, 0] < cx).sum() == 0)
                rows.append({"held": held, "level": lv, "ligament": lig, "n_ref": len(g), "n_inside": int(inside.sum()),
                             "selection_size": len(sel), "projection_size": len(proj),
                             "fallback": int(full.sum() < 64), "side_empty": int(empty)})
    df = pd.DataFrame(rows)
    g = df.groupby("ligament")
    tab = pd.DataFrame({"ref": g.n_ref.sum(), "inside_pct": (g.n_inside.sum() / g.n_ref.sum() * 100).round(1),
                        "selection_median": g.selection_size.median(), "projection_median": g.projection_size.median(),
                        "fallback": g.fallback.sum(), "side_empty": g.side_empty.sum()}).reindex(LIGAMENTS)
    tab.to_csv(OUT / "table_S4_candidate_region_feasibility.csv")
    res["zones_inside_pct"] = r(df.n_inside.sum() / df.n_ref.sum() * 100, 1)
    res["zones_inside_pct_by_class"] = tab.inside_pct.to_dict()
    res["zones_inside_pct_by_level"] = {k: r(v.n_inside.sum() / v.n_ref.sum() * 100, 1) for k, v in df.groupby("level")}


def development(res):
    _, held = read_manifest()
    gt = {k: g for k in held if (g := cases.landmarks(*k)) is not None}
    ck = pickle.loads((P1 / "checkpoints" / "gcz_input_for_gczrl.pkl").read_bytes())
    s = summarize(evaluate(ck["held"], gt, "GCZ input"), LIGAMENTS)
    res["dev_GCZ_input_to_GCZRL"] = {"macro_mae": r(s["macro_mae_mm"]), "pass5_pooled": r(s["pass5"], 1)}
    m = pd.read_csv(P1 / "outputs" / "per_vertebra_per_ligament_metrics.csv")
    rl = m[(m.method == "GCZ-RL") & ~m.missing_prediction].set_index(["case", "ligament"]).mae_mm
    inp = evaluate(ck["held"], gt, "GCZ input")
    inp = inp[~inp.missing_prediction].set_index(["case", "ligament"]).mae_mm
    j = pd.concat([inp.rename("input"), rl.rename("rl")], axis=1).dropna()
    res["dev_GCZRL_vs_its_input"] = {"cases": int(len(j)), "wilcoxon_p": float(wilcoxon(j.rl, j["input"]).pvalue)}


def check(res):
    exp = json.loads((REPO / "paper" / "expected_values.json").read_text())
    fails = []

    def get(d, path):
        for k in path.split("."):
            d = d[k]
        return d
    for path, want in exp.items():
        got = get(res, path)
        if isinstance(want, list):
            ok = len(got) == len(want) and all(abs(a - b) <= 0.0006 for a, b in zip(got, want))
        elif isinstance(want, float):
            ok = abs(got - want) <= 0.0006 + 1e-9
        else:
            ok = got == want
        if not ok:
            fails.append(f"{path}: expected {want}, got {got}")
    for f in fails:
        print("MISMATCH", f)
    print(f"{len(exp) - len(fails)}/{len(exp)} manuscript values reproduced")
    return not fails


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    res = {}
    primary(res)
    zones(res)
    development(res)
    (OUT / "paper_numbers.json").write_text(json.dumps(res, indent=2), encoding="utf-8")
    sys.exit(0 if check(res) else 1)


if __name__ == "__main__":
    main()
