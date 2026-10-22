"""Part 2 - Patient-level validation: leave-one-subject-out (LOSO) on 14 lumbar vertebrae.

Cohort (every annotated lumbar vertebra in the dataset)
    S01: L1-L5 (5)    S02: L1-L4 (4; the S02 scan ends at L4)    S03: L1-L5 (5)

For each of the three folds one whole subject is held out.  Using the two training subjects
only, the script
    1. builds the GCZ zones (thresholds) and uses the training vertebrae as GCZ atlases;
    2. produces training-side GCZ predictions: each training vertebra is predicted with the
       zones and atlases of the *other* training subject, as a test vertebra would be;
    3. builds residual training samples and correction radii;
    4. selects, per ligament class, the Extra Trees configuration (ET-1..ET-4) and the blend
       weight alpha by an inner leave-one-subject-out on the two training subjects;
    5. trains the selected Extra Trees on both training subjects and predicts the held-out
       subject.
The held-out subject's annotations are read only in the final evaluation step.
SLD (official SpineToolkit code) has no training step and is run once per subject.

Outputs (outputs/)
    table_P2_1_folds.csv             fold composition and training sample counts per ligament
    table_P2_2_overall.csv           MAE, Pass@5 and coverage metrics per method
    table_P2_3_per_ligament.csv      MAE and Pass@5 per ligament
    table_P2_4_per_subject.csv       MAE and Pass@5 per held-out subject
    table_P2_6_selected_models.csv   Extra Trees configuration and alpha chosen in each fold
    figure_P2_1_per_ligament.png/pdf per-ligament MAE
    predictions_*.csv, models/, fold_details.json, run_log_part2.txt

Usage:  python run_part2.py            (GCZ geometry and SLD output from checkpoints)
        python run_part2.py --fresh    (recompute them as well, about 25 min extra)
"""
from __future__ import annotations

import argparse
import json
import pickle
import sys
import time
from pathlib import Path

import joblib
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT / "code" / "ligament"))

import cases  # noqa: E402
import gcz  # noqa: E402
import gczrl  # noqa: E402
from data_io import LIGAMENTS, RAW_ROOT, SUBJECT_CODE, case_code, read_vtk_polydata  # noqa: E402
from metrics import evaluate, summarize, to_rows  # noqa: E402

OUT, CKPT = HERE / "outputs", HERE / "checkpoints"
COHORT = {"Verse005": ["L1", "L2", "L3", "L4", "L5"],
          "Verse561": ["L1", "L2", "L3", "L4"],
          "Verse581": ["L1", "L2", "L3", "L4", "L5"]}
SURFACE = "vtk"
ET_CONFIGS = ["ET-1", "ET-2", "ET-3", "ET-4"]
ALPHAS = [0.0, 0.25, 0.40, 0.55, 0.70, 0.85, 1.00, 1.10, 1.20]
gczrl.SIDE_LOCKED = {"ITL", "CL"}   # lateral classes keep the left/right side of their GCZ point
LOG = []


def log(msg):
    line = f"[{time.strftime('%H:%M:%S')}] {msg}"
    print(line, flush=True)
    LOG.append(line)


def keys_of(subjects):
    return [(s, lv) for s in subjects for lv in COHORT[s]]


def gcz_predict(targets, atlas_keys, zone_keys, seed0):
    zone_data = cases.train_dict(zone_keys, SURFACE)
    thr = gcz.build_thresholds(zone_data)
    ref = next((k for k in zone_keys if k[1] == "L3"), zone_keys[0])
    fb = gcz.build_fallback_template(zone_data, ref, dict(n_sub=5000, max_iter=100, tol=1e-6, seed=42))
    atlas = {k: {"verts": cases.mesh(*k, SURFACE), "gt": cases.landmarks(*k)} for k in atlas_keys}
    out, chosen = {}, {}
    for i, key in enumerate(targets):
        s = seed0 + i
        p, best, res = gcz.predict_case(cases.mesh(*key, SURFACE), atlas, thr, fb,
                                        quick_seeds=[300 + s + j for j in range(len(atlas))],
                                        sample_seed=400 + s, full_seed=200 + s, full_tol=1e-6)
        out[key], chosen[key] = p, f"{case_code(*best)} (residual {res:.2f} mm)"
    return out, chosen


def gcz_geometry(test_subject, train_subjects, fresh):
    tag = f"fold_held_{SUBJECT_CODE[test_subject]}"
    ck = CKPT / f"{tag}_gcz.pkl"
    if ck.exists() and not fresh:
        log(f"{tag}: GCZ geometry read from checkpoint (use --fresh to recompute)")
        return pickle.loads(ck.read_bytes()), "checkpoint"
    train_keys, test_keys = keys_of(train_subjects), keys_of([test_subject])
    test_gcz, test_atlas = gcz_predict(test_keys, train_keys, train_keys, seed0=1000)
    train_gcz, train_atlas = {}, {}
    for s in train_subjects:
        other = [k for k in train_keys if k[0] != s]
        p, ch = gcz_predict(keys_of([s]), other, other, seed0=2000 + 100 * len(train_gcz))
        train_gcz.update(p)
        train_atlas.update(ch)
    geo = {"test_gcz": test_gcz, "train_gcz": train_gcz, "test_atlas": test_atlas, "train_atlas": train_atlas}
    ck.write_bytes(pickle.dumps(geo))
    return geo, "computed"


def run_fold(test_subject, fresh):
    train_subjects = [s for s in COHORT if s != test_subject]
    train_keys, test_keys = keys_of(train_subjects), keys_of([test_subject])
    tag = f"fold_held_{SUBJECT_CODE[test_subject]}"
    geo, geo_state = gcz_geometry(test_subject, train_subjects, fresh)

    verts_of = lambda k: cases.mesh(*k, SURFACE)
    labels = lambda k: cases.landmarks(*k) if k in train_keys else None   # held-out labels unreachable
    train = {k: {"verts": verts_of(k), "gt": labels(k)} for k in train_keys}
    ctx = gczrl.make_context(train)
    gcz_train = lambda k: geo["train_gcz"].get(k)
    X, Y, meta = gczrl.training_samples(ctx, train_keys, verts_of, gcz_train, labels)
    radius = gczrl.correction_radius(Y, meta)
    counts = {lig: sum(1 for m in meta if m[2] == lig) for lig in LIGAMENTS}

    inner_ctx, inner_rad = {}, {}
    for s in train_subjects:
        inner_ctx[s] = gczrl.make_context({k: v for k, v in train.items() if k[0] != s})
        rows = [i for i in range(len(meta)) if meta[i][0] != s]
        inner_rad[s] = gczrl.correction_radius(Y[rows], [meta[i] for i in rows])
    grid, models = [], {}
    for name in ET_CONFIGS:
        rec = gczrl.held_out_records(ctx, X, Y, meta, name, radius, train_keys, verts_of, gcz_train, labels,
                                     group_of=lambda k: k[0], ctx_without=inner_ctx.get,
                                     radius_without=inner_rad.get, original_rows_only=True)
        for lig in LIGAMENTS:
            for a in ALPHAS:
                grid.append({"fold": tag, "regressor": name, "ligament": lig, "alpha": a,
                             "inner_loso_mae_mm": gczrl.record_error(ctx, rec[lig], a, lig)})
    grid = pd.DataFrame(grid)
    selection = {}
    for lig in LIGAMENTS:
        g = grid[grid.ligament == lig].sort_values(["inner_loso_mae_mm", "regressor", "alpha"]).iloc[0]
        selection[lig] = (g.regressor, float(g.alpha), float(g.inner_loso_mae_mm))
    for name in sorted({v[0] for v in selection.values()}):
        models[name] = gczrl.fit(name, X, Y)
    pred = {}
    for key in test_keys:
        verts, base = verts_of(key), geo["test_gcz"][key] or {}
        pred[key] = {}
        for lig in LIGAMENTS:
            gp = base.get(lig)
            if gp is None:
                continue
            name, alpha, _ = selection[lig]
            corr = gczrl.corrected_points(ctx, models[name], radius[lig], key, verts, gp, lig)
            pred[key][lig] = gczrl.blend(ctx, gp, corr, alpha, verts, lig)
    joblib.dump({"models": models, "selection": selection, "correction_radius_mm": radius,
                 "thresholds": ctx.thresholds, "benchmark_prior": ctx.bench_norm},
                OUT / "models" / f"{tag}_extra_trees.joblib", compress=3)
    log(f"{tag}: {len(train_keys)} training vertebrae, {len(X)} samples; selected "
        + ", ".join(f"{l}={selection[l][0]}/a={selection[l][1]}" for l in LIGAMENTS))
    info = {"fold": tag, "held_out_subject": SUBJECT_CODE[test_subject],
            "training_subjects": [SUBJECT_CODE[s] for s in train_subjects],
            "training_vertebrae": [case_code(*k) for k in train_keys],
            "test_vertebrae": [case_code(*k) for k in test_keys],
            "training_samples": int(len(X)), "training_samples_per_ligament": counts,
            "selection": {l: {"extra_trees": v[0], "alpha": v[1], "inner_loso_mae_mm": v[2]} for l, v in selection.items()},
            "correction_radius_mm": radius, "gcz_geometry": geo_state,
            "test_atlas_choice": {case_code(*k): v for k, v in geo["test_atlas"].items()},
            "held_out_annotations_read_before_evaluation": False}
    return geo["test_gcz"], pred, info, grid


def run_sld(fresh):
    ck = CKPT / "sld_official_medial.pkl"
    if ck.exists() and not fresh:
        log("SLD output read from checkpoint (use --fresh to recompute)")
        return pickle.loads(ck.read_bytes()), "checkpoint"
    import sld_official
    preds, diag = {}, {}
    for subj, levels in COHORT.items():
        out = sld_official.run_subject({lv: read_vtk_polydata(RAW_ROOT / subj / f"{subj}_{lv}.vtk") for lv in levels},
                                       medial=True)
        d = out.pop("_diag")
        diag[SUBJECT_CODE[subj]] = {k: d[k] for k in ("errors", "seconds", "levels") if k in d}
        for lv in levels:
            preds[(subj, lv)] = out.get(lv, {})
        log(f"SLD {SUBJECT_CODE[subj]}: {d['seconds']} s, errors: {len(d['errors'])}")
    res = {"pred": preds, "diag": diag}
    ck.write_bytes(pickle.dumps(res))
    return res, "computed"


def figure(per_lig):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    cols = {"GCZ": "#2a78d6", "GCZ-RL": "#eb6834", "SLD": "#1baf7a"}
    x = np.arange(len(LIGAMENTS))
    fig, ax = plt.subplots(figsize=(7, 2.8))
    for i, m in enumerate(cols):
        vals = [float(per_lig[(per_lig.method == m) & (per_lig.ligament == l)].MAE_mm.iloc[0]) for l in LIGAMENTS]
        ax.bar(x + (i - 1) * 0.27, vals, width=0.25, color=cols[m], label=m, edgecolor="white", linewidth=0.5)
    ax.axhline(5, color="#52514e", linestyle=(0, (3, 2)), linewidth=0.8)
    ax.set_xticks(x, LIGAMENTS)
    ax.set_ylabel("MAE (mm)")
    ax.spines[["top", "right"]].set_visible(False)
    ax.yaxis.grid(True, color="#e4e3df", linewidth=0.6)
    ax.set_axisbelow(True)
    ax.legend(ncol=3, frameon=False, loc="upper left")
    fig.tight_layout()
    fig.savefig(OUT / "figure_P2_1_per_ligament.png", dpi=200)
    fig.savefig(OUT / "figure_P2_1_per_ligament.pdf")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fresh", action="store_true")
    args = ap.parse_args()
    for d in (OUT, CKPT, OUT / "models"):
        d.mkdir(exist_ok=True)
    t0 = time.time()
    pred_gcz, pred_rl, folds, grids = {}, {}, [], []
    for subj in COHORT:
        g, r, info, grid = run_fold(subj, args.fresh)
        pred_gcz.update(g)
        pred_rl.update(r)
        folds.append(info)
        grids.append(grid)
    sld, sld_state = run_sld(args.fresh)

    # ---------- evaluation (the only place where held-out annotations are read)
    all_keys = keys_of(COHORT)
    gt = {k: cases.landmarks(*k) for k in all_keys}
    code = lambda k: SUBJECT_CODE[k[0]]
    methods = {"GCZ": pred_gcz, "GCZ-RL": pred_rl, "SLD": sld["pred"]}
    overall, lig_rows, subj_rows, cl = [], [], [], []
    for m, pred in methods.items():
        df = evaluate(pred, gt, m, subject_of=code)
        cl.append(df)
        s = summarize(df, LIGAMENTS)
        overall.append({"method": m, "vertebrae": s["cases"], "MAE_mm": s["macro_mae_mm"], "Pass5_percent": s["pass5"],
                        "reference_coverage_percent": s["ref_coverage_mean"] * 100,
                        "symmetric_nn_distance_mm": s["symmetric_nn_mm"],
                        "all_reference_points_within_5mm_percent": s["ref_pass5_all_mean"],
                        "missing_ligament_outputs": s["missing_case_ligaments"]})
        for lig in LIGAMENTS:
            lig_rows.append({"method": m, "ligament": lig, "MAE_mm": s["per_ligament_mae_mm"][lig],
                             "Pass5_percent": s["per_ligament_pass5"][lig]})
        for subj, sub in df.groupby("subject"):
            ss = summarize(sub, LIGAMENTS)
            subj_rows.append({"method": m, "held_out_subject": subj, "vertebrae": ss["cases"],
                              "MAE_mm": ss["macro_mae_mm"], "Pass5_percent": ss["pass5"]})
        to_rows(pred, m, code=lambda p, l: case_code(p, l)).to_csv(
            OUT / f"predictions_{m.replace('-', '').lower()}.csv", index=False)
    pd.DataFrame(overall).round(3).to_csv(OUT / "table_P2_2_overall.csv", index=False)
    per_lig = pd.DataFrame(lig_rows)
    per_lig.round(3).to_csv(OUT / "table_P2_3_per_ligament.csv", index=False)
    pd.DataFrame(subj_rows).round(3).to_csv(OUT / "table_P2_4_per_subject.csv", index=False)
    pd.concat(cl).to_csv(OUT / "per_vertebra_per_ligament_metrics.csv", index=False)
    pd.DataFrame([{"held_out_subject": f["held_out_subject"], "training_subjects": " + ".join(f["training_subjects"]),
                   "training_vertebrae": len(f["training_vertebrae"]), "test_vertebrae": len(f["test_vertebrae"]),
                   "training_samples": f["training_samples"], **{f"samples_{l}": f["training_samples_per_ligament"][l] for l in LIGAMENTS}}
                  for f in folds]).to_csv(OUT / "table_P2_1_folds.csv", index=False)
    pd.DataFrame([{"held_out_subject": f["held_out_subject"], "ligament": l, **f["selection"][l]}
                  for f in folds for l in LIGAMENTS]).round(3).to_csv(OUT / "table_P2_6_selected_models.csv", index=False)
    pd.concat(grids).to_csv(OUT / "inner_selection_grid.csv", index=False)
    figure(per_lig)
    json.dump(folds, open(OUT / "fold_details.json", "w"), indent=2)
    json.dump({"sld_state": sld_state, "sld_diagnostics": sld["diag"], "folds_completed": len(folds),
               "seconds": round(time.time() - t0)}, open(OUT / "run_status_part2.json", "w"), indent=2)
    (OUT / "run_log_part2.txt").write_text("\n".join(LOG), encoding="utf-8")
    log("done: " + "; ".join(f"{r['method']} {r['MAE_mm']:.3f} mm / {r['Pass5_percent']:.1f}%" for r in overall))


if __name__ == "__main__":
    main()
