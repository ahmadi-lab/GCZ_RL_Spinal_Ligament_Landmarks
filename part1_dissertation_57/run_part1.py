"""Part 1 - Dissertation Chapter 4: 57 vertebrae, 20 calibration / 37 held-out.

Produces every data table/figure of Chapter 4 from the raw meshes and annotations:

  Table 4.1   dataset composition                -> outputs/table_4_1_dataset.csv
  Table 4.3   calibration comparison of regressors -> outputs/table_4_3_regressor_comparison.csv
  Table 4.4   held-out MAE / Pass@5 per ligament  -> outputs/table_4_4_heldout_results.csv
  Figure 4.10 bar chart of Table 4.4              -> outputs/figure_4_10_heldout_results.png/.pdf
  Figure 4.11 landmarks on a held-out vertebra    -> outputs/figure_4_11_landmarks_*.png + .csv
  check file  regenerated vs dissertation values  -> outputs/verification_part1.json

Steps
  1. GCZ on the 37 held-out vertebrae (atlas = the 5 S01 lumbar calibration vertebrae).
  2. GCZ on the 20 calibration vertebrae (an atlas vertebra never predicts itself) and a
     second GCZ run on the held-out vertebrae that serves as GCZ-RL input.
  3. Residual training samples from the 20 calibration vertebrae (1683 samples x 36 features).
  4. Table 4.3: leave-one-vertebra-out on the 20 calibration vertebrae for every candidate
     regressor and blend weight.
  5. Final GCZ-RL (Extra Trees, configuration in code/ligament/gczrl.py) trained on all 20
     calibration vertebrae and applied once to the held-out vertebrae.
  6. Evaluation against the expert annotations of the 27 held-out VerSe vertebrae, as in the
     dissertation (all 57 vertebrae are annotated; the 10 cervical ones were not scored).

Usage:  python run_part1.py            (GCZ geometry from checkpoints, everything else recomputed)
        python run_part1.py --fresh    (recompute the GCZ geometry as well, about 40 min extra)
"""
from __future__ import annotations

import argparse
import hashlib
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
from data_io import DATASET, LIGAMENTS, RAW_BASE, SUBJECT_CODE, case_code, landmark_dir, read_manifest  # noqa: E402
from metrics import coordinate_agreement, evaluate, summarize, table_from_csv, to_rows  # noqa: E402

OUT, CKPT = HERE / "outputs", HERE / "checkpoints"
ARCHIVED = DATASET / "dissertation_archived_predictions"
ATLAS = [("Verse005", l) for l in ["L1", "L2", "L3", "L4", "L5"]]
ALPHAS = [0.70, 0.85, 1.00, 1.10, 1.20]
TABLE_4_3_MODELS = {"Ridge": ["Ridge"], "KNN": ["KNN"], "Random forest": ["RandomForest"],
                    "Gradient boosting": ["GradientBoosting", "HistGradientBoosting"],
                    "Extra trees": ["ET-1", "ET-2", "ET-3", "ET-4"]}
DISSERTATION_TABLE_4_4 = {  # values printed in the dissertation, used only for the check file
    "GCZ": {"ALL": 3.050, "PLL": 2.329, "LF": 4.110, "ISL": 2.495, "SSL": 2.921, "ITL": 5.614, "CL": 5.437, "Mean": 3.708},
    "GCZ-RL": {"ALL": 2.708, "PLL": 2.164, "LF": 2.933, "ISL": 2.277, "SSL": 2.681, "ITL": 4.434, "CL": 3.426, "Mean": 2.946},
    "SLD": {"ALL": 2.750, "PLL": 2.246, "LF": 4.814, "ISL": 4.560, "SSL": 1.739, "ITL": 3.485, "CL": 1.630, "Mean": 3.032}}
DISSERTATION_PASS = {
    "GCZ": {"ALL": 86.7, "PLL": 92.6, "LF": 67.0, "ISL": 91.2, "SSL": 85.2, "ITL": 72.2, "CL": 63.7, "Mean": 79.8},
    "GCZ-RL": {"ALL": 88.1, "PLL": 94.4, "LF": 88.5, "ISL": 94.6, "SSL": 88.9, "ITL": 81.5, "CL": 80.7, "Mean": 88.1},
    "SLD": {"ALL": 90.4, "PLL": 91.8, "LF": 60.9, "ISL": 69.5, "SSL": 100.0, "ITL": 76.7, "CL": 99.0, "Mean": 84.0}}
LOG = []


def log(msg):
    line = f"[{time.strftime('%H:%M:%S')}] {msg}"
    print(line, flush=True)
    LOG.append(line)


def cached(name, fn, fresh):
    path = CKPT / f"{name}.pkl"
    if path.exists() and not fresh:
        log(f"GCZ geometry read from checkpoint {path.name} (use --fresh to recompute)")
        return pickle.loads(path.read_bytes())
    t0 = time.time()
    obj = fn()
    path.write_bytes(pickle.dumps(obj))
    log(f"GCZ geometry computed in {time.time() - t0:.0f} s -> {path.name}")
    return obj


# ---------------------------------------------------------------------------- steps 1-2
def gcz_heldout_for_table(calib, held):
    """Step 1: GCZ column of Table 4.4 (STL surfaces, seeds as in the original run)."""
    cd = cases.train_dict(calib, "stl")
    thr = gcz.build_thresholds(cd)
    fb = gcz.build_fallback_template(cd, ("Verse005", "L1"), dict(n_sub=5000, max_iter=100, tol=1e-6, seed=42))
    atlas = {k: cd[k] for k in ATLAS}
    pred = {}
    for k, key in enumerate(held):
        pred[key], _, _ = gcz.predict_case(cases.mesh(*key, "stl"), atlas, thr, fb,
                                           quick_seeds=[300 + k + a for a in range(5)],
                                           sample_seed=400 + k, full_seed=200 + k, full_tol=1e-6)
    return {"pred": pred, "thresholds": thr}


def gcz_input_for_gczrl(calib, held):
    """Step 2: GCZ predictions that GCZ-RL corrects (VTK surfaces).

    Calibration vertebrae: an S01 atlas vertebra is predicted from the other four atlas
    vertebrae, so no vertebra is predicted from its own annotations."""
    cd = cases.train_dict(calib, "vtk")
    thr = gcz.build_thresholds(cd)
    ref = next(c for c in calib if c[1] == "L3")
    fb = gcz.build_fallback_template(cd, ref, dict(n_sub=5000, max_iter=100, seed=42))

    def one(key, atlas_keys, seed):
        atlas = {k: cd[k] for k in atlas_keys}
        p, _, _ = gcz.predict_case(cases.mesh(*key, "vtk"), atlas, thr, fb,
                                   quick_seeds=[300 + j + seed for j in range(len(atlas))],
                                   sample_seed=seed, full_seed=200 + seed, full_tol=1e-5, third_fallback=False)
        return p
    calib_pred = {key: one(key, [a for a in ATLAS if a != key], 42 + i * 7) for i, key in enumerate(calib)}
    held_pred = {key: one(key, ATLAS, 42 + k) for k, key in enumerate(held)}
    return {"calib": calib_pred, "held": held_pred, "thresholds": thr}


# ------------------------------------------------------------------------------ helpers
def code_of(key):
    return SUBJECT_CODE.get(key[0], key[0])


def per_ligament(pred, gt, name):
    df = evaluate(pred, gt, name, subject_of=code_of)
    s = summarize(df, LIGAMENTS)
    return df, s


CERVICAL_GT = RAW_BASE / "Cervical_Ligaments_GT"   # public SpineToolkit cervical annotations


def has_annotation(key):
    """Expert annotations exist for every vertebra: VerSe folders for S01-S03, the
    Cervical_Ligaments_GT folder for S04-S05."""
    p, lv = key
    return landmark_dir(p, lv) is not None or (CERVICAL_GT / p / f"{p}_{lv}_ligaments").is_dir()


def table_4_1(calib, held):
    """Scored = annotations used by the original evaluation (VerSe folders only); the
    cervical annotations were not loaded by the original notebook (BASE_CERV = BASE)."""
    rows = []
    for split, keys in (("Calibration", calib), ("Held-out test", held)):
        ann = [k for k in keys if has_annotation(k)]
        scored = [k for k in keys if landmark_dir(*k) is not None]
        rows.append({"split": split, "vertebrae": len(keys), "with_expert_annotations": len(ann),
                     "annotations_used_in_evaluation": len(scored),
                     "subjects": ", ".join(sorted({code_of(k) for k in keys}))})
    rows.append({"split": "Total", "vertebrae": len(calib) + len(held),
                 "with_expert_annotations": sum(r["with_expert_annotations"] for r in rows),
                 "annotations_used_in_evaluation": sum(r["annotations_used_in_evaluation"] for r in rows),
                 "subjects": "S01-S05"})
    detail = pd.DataFrame([{"case": case_code(*k), "subject": code_of(k), "level": k[1],
                            "split": "calibration" if k in calib else "held-out",
                            "has_annotations": has_annotation(k),
                            "annotations_used_in_evaluation": landmark_dir(*k) is not None} for k in calib + held])
    return pd.DataFrame(rows), detail


def figure_4_10(t44):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    cols = {"GCZ": "#2a78d6", "GCZ-RL": "#eb6834", "SLD": "#1baf7a"}
    labels = LIGAMENTS + ["Mean"]
    x = np.arange(len(labels))
    fig, axes = plt.subplots(2, 1, figsize=(7.5, 5.2), sharex=True)
    for ax, metric, ylab in [(axes[0], "MAE_mm", "MAE (mm)"), (axes[1], "Pass5_percent", "Pass@5 mm (%)")]:
        for i, m in enumerate(["GCZ", "GCZ-RL", "SLD"]):
            vals = [float(t44[(t44.method == m) & (t44.ligament == l)][metric].iloc[0]) for l in labels]
            ax.bar(x + (i - 1) * 0.27, vals, width=0.25, color=cols[m], label=m, edgecolor="white", linewidth=0.5)
        ax.set_ylabel(ylab)
        ax.spines[["top", "right"]].set_visible(False)
        ax.yaxis.grid(True, color="#e4e3df", linewidth=0.6)
        ax.set_axisbelow(True)
    axes[0].axhline(5, color="#52514e", linestyle=(0, (3, 2)), linewidth=0.8)
    axes[0].legend(ncol=3, frameon=False, loc="upper left")
    axes[1].set_xticks(x, labels)
    fig.tight_layout()
    fig.savefig(OUT / "figure_4_10_heldout_results.png", dpi=200)
    fig.savefig(OUT / "figure_4_10_heldout_results.pdf")
    plt.close(fig)


def figure_4_11(pred, key):
    """Landmarks predicted by GCZ-RL on one held-out vertebra, coloured by ligament class."""
    import pyvista as pv
    from data_io import RAW_ROOT, read_vtk_polydata
    cols = {"ALL": "red", "PLL": "limegreen", "LF": "royalblue", "ISL": "orange",
            "SSL": "purple", "ITL": "magenta", "CL": "deepskyblue"}
    v, f = read_vtk_polydata(RAW_ROOT / key[0] / f"{key[0]}_{key[1]}.vtk")
    mesh = pv.PolyData(v, np.hstack([np.full((len(f), 1), 3), f]).ravel())
    pl = pv.Plotter(off_screen=True, window_size=(900, 800))
    pl.set_background("white")
    pl.add_mesh(mesh, color="wheat", opacity=0.55)
    for lig, pts in pred[key].items():
        pl.add_points(np.atleast_2d(pts), color=cols[lig], point_size=14, render_points_as_spheres=True, label=lig)
    pl.add_legend(bcolor="white", face="circle", size=(0.18, 0.28))
    pl.view_vector((0, 0, 1), viewup=(0, -1, 0))   # axial view from superior, posterior up
    name = f"figure_4_11_landmarks_{case_code(*key)}.png"
    pl.screenshot(OUT / name)
    pl.close()
    to_rows({key: pred[key]}, "GCZ-RL", code=lambda p, l: case_code(p, l)).to_csv(
        OUT / f"figure_4_11_landmarks_{case_code(*key)}.csv", index=False)
    return name


# --------------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fresh", action="store_true", help="recompute GCZ geometry")
    args = ap.parse_args()
    OUT.mkdir(exist_ok=True)
    CKPT.mkdir(exist_ok=True)
    (OUT / "models").mkdir(exist_ok=True)
    t0 = time.time()
    calib, held = read_manifest()
    gt_held = {k: g for k in held if (g := cases.landmarks(*k)) is not None}
    log(f"57 vertebrae: {len(calib)} calibration, {len(held)} held-out; "
        f"{len(gt_held)} held-out vertebrae have expert annotations")

    t41, t41_detail = table_4_1(calib, held)
    t41.to_csv(OUT / "table_4_1_dataset.csv", index=False)
    t41_detail.to_csv(OUT / "table_4_1_dataset_per_vertebra.csv", index=False)

    step1 = cached("gcz_heldout_table44", lambda: gcz_heldout_for_table(calib, held), args.fresh)
    step2 = cached("gcz_input_for_gczrl", lambda: gcz_input_for_gczrl(calib, held), args.fresh)

    # ---- step 3: training samples (calibration vertebrae only)
    train = {k: {"verts": cases.mesh(*k, "vtk"), "gt": cases.landmarks(*k)} for k in calib}
    ctx = gczrl.make_context(train)
    verts_of = lambda k: cases.mesh(*k, "vtk")
    gt_calib = lambda k: cases.landmarks(*k) if k in calib else None   # held-out labels unreachable
    X, Y, meta = gczrl.training_samples(ctx, calib, verts_of, lambda k: step2["calib"].get(k), gt_calib)
    radius = gczrl.correction_radius(Y, meta)
    log(f"training samples: {len(X)} x {X.shape[1]} from {len(calib)} calibration vertebrae")

    # ---- step 4: Table 4.3 (leave-one-vertebra-out on calibration vertebrae)
    grid = []
    for family, names in TABLE_4_3_MODELS.items():
        for name in names:
            rec = gczrl.held_out_records(ctx, X, Y, meta, name, radius, calib, verts_of,
                                         lambda k: step2["calib"].get(k), gt_calib)
            for lig in LIGAMENTS:
                for a in ALPHAS:
                    grid.append({"family": family, "regressor": name, "ligament": lig, "alpha": a,
                                 "calibration_mae_mm": gczrl.record_error(ctx, rec[lig], a, lig),
                                 "vertebrae": len(rec[lig])})
            log(f"Table 4.3: {name} evaluated (leave-one-vertebra-out, 20 folds)")
    grid = pd.DataFrame(grid)
    grid.to_csv(OUT / "table_4_3_full_grid.csv", index=False)
    best = grid.groupby(["ligament", "family"]).calibration_mae_mm.min().unstack()[list(TABLE_4_3_MODELS)]
    best = best.reindex(LIGAMENTS)
    best["lowest"] = best.idxmin(axis=1)
    best.round(3).to_csv(OUT / "table_4_3_regressor_comparison.csv")

    # ---- step 5: final GCZ-RL (Extra Trees)
    final_models, final = {}, {}
    for name in sorted({m for m, _ in gczrl.FINAL_CONFIGURATION.values()}):
        final_models[name] = gczrl.fit(name, X, Y)
    joblib.dump({"models": final_models, "configuration": gczrl.FINAL_CONFIGURATION, "correction_radius_mm": radius,
                 "thresholds": ctx.thresholds, "benchmark_prior": ctx.bench_norm},
                OUT / "models" / "gczrl_extra_trees_final.joblib", compress=3)
    for key in held:
        verts, base = verts_of(key), step2["held"][key] or {}
        final[key] = {}
        for lig in LIGAMENTS:
            gp = base.get(lig)
            if gp is None:
                continue
            name, alpha = gczrl.FINAL_CONFIGURATION[lig]
            corr = gczrl.corrected_points(ctx, final_models[name], radius[lig], key, verts, gp, lig)
            final[key][lig] = gczrl.blend(ctx, gp, corr, alpha, verts, lig)
    log("final GCZ-RL (Extra Trees) trained on 20 calibration vertebrae and applied to 37 held-out vertebrae")

    # ---- step 6: evaluation
    sld = table_from_csv(ARCHIVED / "sld_heldout37.csv")
    methods = {"GCZ": step1["pred"], "GCZ-RL": final, "SLD": sld}
    rows, cl = [], []
    for m, pred in methods.items():
        df, s = per_ligament(pred, gt_held, m)
        cl.append(df)
        for lig in LIGAMENTS:
            rows.append({"method": m, "ligament": lig, "MAE_mm": s["per_ligament_mae_mm"][lig],
                         "Pass5_percent": s["per_ligament_pass5"][lig]})
        rows.append({"method": m, "ligament": "Mean", "MAE_mm": s["macro_mae_mm"], "Pass5_percent": s["pass5"]})
    t44 = pd.DataFrame(rows)
    t44.round(3).to_csv(OUT / "table_4_4_heldout_results.csv", index=False)
    pd.concat(cl).to_csv(OUT / "per_vertebra_per_ligament_metrics.csv", index=False)
    for m, pred in [("gcz", step1["pred"]), ("gczrl_extra_trees", final)]:
        to_rows(pred, m, code=lambda p, l: case_code(p, l)).to_csv(OUT / f"predictions_{m}_heldout37.csv", index=False)
    figure_4_10(t44)
    fig411 = figure_4_11(final, ("Verse581", "L3"))

    # ---- checks
    agree = {m: coordinate_agreement(pred, table_from_csv(ARCHIVED / f))
             for m, pred, f in [("GCZ", step1["pred"], "gcz_heldout37.csv"),
                                ("GCZ-RL", final, "gczrl_extra_trees_heldout37.csv")]}
    diffs = {m: {l: round(float(t44[(t44.method == m) & (t44.ligament == l)].MAE_mm.iloc[0]), 3) - v
                 for l, v in DISSERTATION_TABLE_4_4[m].items()} for m in DISSERTATION_TABLE_4_4}
    pdiffs = {m: {l: round(float(t44[(t44.method == m) & (t44.ligament == l)].Pass5_percent.iloc[0]), 1) - v
                  for l, v in DISSERTATION_PASS[m].items()} for m in DISSERTATION_PASS}
    check = {
        "regenerated_vs_archived_predictions": agree,
        "table_4_4_MAE_minus_dissertation": diffs,
        "table_4_4_Pass5_minus_dissertation": pdiffs,
        "table_4_4_matches_dissertation": all(abs(d) < 1e-9 for m in diffs.values() for d in m.values())
        and all(abs(d) < 1e-9 for m in pdiffs.values() for d in m.values()),
        "training_vertebrae": [case_code(*k) for k in calib],
        "held_out_annotations_used_for_training": False,
        "training_samples": int(len(X)), "features": int(X.shape[1]),
        "figure_4_11_file": fig411, "seconds": round(time.time() - t0),
        "sha256_final_model": hashlib.sha256((OUT / "models" / "gczrl_extra_trees_final.joblib").read_bytes()).hexdigest(),
    }
    (OUT / "verification_part1.json").write_text(json.dumps(check, indent=2), encoding="utf-8")
    (OUT / "run_log_part1.txt").write_text("\n".join(LOG), encoding="utf-8")
    log(f"Table 4.4 reproduced exactly: {check['table_4_4_matches_dissertation']}; "
        f"GCZ-RL identical point sets: {agree['GCZ-RL']['identical_sets_tol_1um']}/{agree['GCZ-RL']['case_ligament_sets_compared']}")


if __name__ == "__main__":
    main()
