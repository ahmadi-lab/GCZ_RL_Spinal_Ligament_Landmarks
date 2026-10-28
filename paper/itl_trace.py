"""Step-by-step trace of the GCZ-RL correction of the two ITL landmarks of held-out vertebra S03 L4
(primary evaluation, fold holding out S03), from the saved outputs; nothing is retrained.

The trace repeats the operations of part2_patient_loso_14/run_part2.py (run_fold) and
code/ligament/gczrl.py with the saved fold model:
  1. GCZ landmark (atlas transfer and snapping; predictions_gcz.csv)
  2. raw residual predicted by the saved Extra Trees model of the fold
  3. residual clipped to the correction radius of ITL
  4. corrected point projected onto the ITL projection region, side kept (corrected_points)
  5. blended with the GCZ landmark by the inner-validated weight alpha
  6. projected again onto the region, side kept (blend) = the GCZ-RL landmark
It also projects the GCZ landmark directly ("projection only", no residual), which separates the
part of the change that the geometry alone produces.

Inputs of the computation: the target mesh, the saved GCZ landmarks, and the saved fold model with
its thresholds, class priors, correction radii, and alpha, all derived from the training subjects
S01 and S02. The expert annotations of S03 are read only afterwards, to score each step. The final
positions must equal the saved GCZ-RL landmarks, and the scored errors the metrics file.

Usage: python paper/itl_trace.py
Outputs: paper/outputs/revision/itl_trace_S03_L4.{csv,json}, tab_itl_trace.tex
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment
from scipy.spatial import cKDTree

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "code" / "ligament"))
import cases  # noqa: E402
import gczrl  # noqa: E402
from metrics import table_from_csv  # noqa: E402

P2 = REPO / "part2_patient_loso_14" / "outputs"
OUT = REPO / "paper" / "outputs" / "revision"
KEY, LIG = ("Verse581", "L4"), "ITL"
gczrl.SIDE_LOCKED = {"ITL", "CL"}            # as in run_part2.py


def main():
    saved = joblib.load(P2 / "models" / "fold_held_S03_extra_trees.joblib")
    ctx = gczrl.Context(thresholds=saved["thresholds"], bench_norm=saved["benchmark_prior"])
    name, alpha, _ = saved["selection"][LIG]
    model, radius = saved["models"][name], saved["correction_radius_mm"][LIG]
    verts = np.asarray(cases.mesh(*KEY, "vtk"), float)
    gp = np.atleast_2d(table_from_csv(P2 / "predictions_gcz.csv")[KEY][LIG]).astype(float)
    final_saved = np.atleast_2d(table_from_csv(P2 / "predictions_gczrl.csv")[KEY][LIG]).astype(float)

    # ---- the correction, without any annotation of S03
    feats = np.array([gczrl.point_features(ctx, verts, p, LIG, KEY[1], key=KEY) for p in gp])
    raw = np.asarray(model.predict(feats), float)
    clipped = np.array([gczrl.clip(r, radius) for r in raw])
    corr = gp + clipped
    corr_proj = gczrl.project_to_zone(ctx, verts, corr, LIG, ref=gp)
    blended = gp + alpha * (corr_proj - gp)
    final = gczrl.project_to_zone(ctx, verts, blended, LIG, ref=gp)
    proj_only = gczrl.project_to_zone(ctx, verts, gp, LIG, ref=gp)
    dev = float(np.abs(final - final_saved).max())
    assert dev < 1e-9, dev                                   # the trace reproduces the saved landmarks
    assert np.allclose(corr_proj, gczrl.corrected_points(ctx, model, radius, KEY, verts, gp, LIG))

    # ---- scoring with the expert annotations of S03
    ref = np.atleast_2d(cases.landmarks(*KEY)[LIG]).astype(float)
    ri, ci = linear_sum_assignment(np.linalg.norm(final[:, None] - ref[None], axis=2))
    pair = dict(zip(ri, ci))
    xc = float(verts[:, 0].mean())                           # midline used by the side constraint
    zone = ctx.zone(verts, LIG)
    in_zone = lambda p: bool(cKDTree(zone).query(p)[0] < 1e-6)  # noqa: E731
    stages = [("1 GCZ landmark", gp), ("2 GCZ + raw residual", gp + raw), ("3 GCZ + clipped residual", corr),
              ("4 projected (corrected point)", corr_proj), (f"5 blended (alpha = {alpha:g})", blended),
              ("6 projected = GCZ-RL", final), ("projection only (no residual)", proj_only)]
    rows = []
    for i in range(len(gp)):
        r = ref[pair[i]]
        side = "left" if gp[i, 0] >= xc else "right"         # LPS: +x is the patient's left
        prev = gp[i]
        for stage, P in stages:
            p = P[i]
            rows.append({"side": side, "stage": stage, "x": p[0], "y": p[1], "z": p[2],
                         "error_mm": float(np.linalg.norm(p - r)),
                         "step_mm": float(np.linalg.norm(p - prev)) if not stage.startswith("projection only") else
                         float(np.linalg.norm(p - gp[i])),
                         "same_side_as_gcz": bool((p[0] >= xc) == (gp[i, 0] >= xc)),
                         "vertex_of_projection_region": in_zone(p)})
            if not stage.startswith("projection only"):
                prev = p
    tab = pd.DataFrame(rows)
    OUT.mkdir(parents=True, exist_ok=True)
    tab.to_csv(OUT / "itl_trace_S03_L4.csv", index=False)                  # full precision (checked against the text)
    met = pd.read_csv(P2 / "per_vertebra_per_ligament_metrics.csv")
    met = met[(met.case == f"{KEY[0]}_{KEY[1]}") & (met.ligament == LIG)].set_index("method").mae_mm
    err = lambda st: float(tab[tab.stage == st].error_mm.mean())  # noqa: E731
    assert abs(err("1 GCZ landmark") - met["GCZ"]) < 1e-9 and abs(err("6 projected = GCZ-RL") - met["GCZ-RL"]) < 1e-9
    summ = {"vertebra": "S03 L4 (Verse581 L4), fold holding out S03", "class": LIG,
            "regressor": name, "alpha": alpha, "correction_radius_mm": round(radius, 3),
            "final_equals_saved_gczrl_max_deviation_mm": dev,
            "class_mae_mm": {st: round(err(st), 3) for st, _ in stages},
            "class_mae_matches_metrics_file": True,
            "raw_residual_length_mm": {s: round(float(np.linalg.norm(raw[i])), 3)
                                       for i, s in enumerate(tab.side.unique())},
            "inputs_without_S03_annotations": "mesh, saved GCZ landmarks, saved fold model (S01, S02)"}
    (OUT / "itl_trace_S03_L4.json").write_text(json.dumps(summ, indent=2))
    lines = []
    for side in ("right", "left"):
        t = tab[tab.side == side]
        for k, r in enumerate(t.itertuples()):
            stage = (r.stage[2:] if r.stage[0].isdigit() else r.stage).replace(f"alpha = {alpha:g}", f"$\\alpha_L={alpha:g}$")
            lines.append(f"{side if k == 0 else ''} & {stage} & "
                         f"{r.step_mm:.2f} & {r.error_mm:.2f} & {'yes' if r.vertex_of_projection_region else 'no'} \\\\")
        lines.append("\\addlinespace")
    (OUT / "tab_itl_trace.tex").write_text("\n".join(lines[:-1]) + "\n", encoding="utf-8")
    print(tab.round(2).to_string(index=False))
    print(json.dumps(summ, indent=1))


if __name__ == "__main__":
    main()
