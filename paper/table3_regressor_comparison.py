"""Recompute the calibration regressor comparison (paper Table 3, Supplementary Table S6; dissertation Table 4.3).

Variants (all leave-one-vertebra-out over the 20 calibration vertebrae):
  as_published   : package code as run for the dissertation (original_rows_only=False,
                   zones/prior/radius from all 20 calibration vertebrae)
  rows_fixed     : only original (non-mirrored) rows of the held-out vertebra are evaluated,
                   each paired with its own GCZ point
  rows_fixed_ctx : as rows_fixed, and zones, class prior and correction radius are rebuilt
                   without the held-out vertebra (as in the patient-level analysis)
The GCZ inputs are the archived checkpoint (gcz_input_for_gczrl.pkl) in every variant.

Usage: python paper/table3_regressor_comparison.py [out_dir]     (about 25 min)
Paper Table 3 = table3_rows_fixed_ctx.csv; Supplementary Table S6 compares all three variants.
"""
import json
import pickle
import sys
import time
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "paper/outputs"
OUT.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(ROOT / "code" / "ligament"))
import cases  # noqa: E402
import gczrl  # noqa: E402
from data_io import LIGAMENTS, read_manifest  # noqa: E402

ALPHAS = [0.70, 0.85, 1.00, 1.10, 1.20]
FAMILIES = {"Ridge": ["Ridge"], "KNN": ["KNN"], "Random forest": ["RandomForest"],
            "Gradient boosting": ["GradientBoosting", "HistGradientBoosting"],
            "Extra trees": ["ET-1", "ET-2", "ET-3", "ET-4"]}

calib, _ = read_manifest()
step2 = pickle.loads((ROOT / "part1_dissertation_57/checkpoints/gcz_input_for_gczrl.pkl").read_bytes())
verts_of = lambda k: cases.mesh(*k, "vtk")
gt_calib = lambda k: cases.landmarks(*k) if k in calib else None
gcz_of = lambda k: step2["calib"].get(k)
train = {k: {"verts": verts_of(k), "gt": cases.landmarks(*k)} for k in calib}
ctx = gczrl.make_context(train)
X, Y, meta = gczrl.training_samples(ctx, calib, verts_of, gcz_of, gt_calib)
radius = gczrl.correction_radius(Y, meta)
print(f"samples {X.shape}", flush=True)

ctx_wo, rad_wo = {}, {}
for k in calib:
    ctx_wo[k] = gczrl.make_context({kk: v for kk, v in train.items() if kk != k})
    rows = [i for i in range(len(meta)) if (meta[i][0], meta[i][1]) != k]
    rad_wo[k] = gczrl.correction_radius(Y[rows], [meta[i] for i in rows])

variants = {
    "as_published": dict(original_rows_only=False),
    "rows_fixed": dict(original_rows_only=True),
    "rows_fixed_ctx": dict(original_rows_only=True, ctx_without=ctx_wo.get, radius_without=rad_wo.get),
}
grid = []
for vname, kw in variants.items():
    for fam, names in FAMILIES.items():
        for name in names:
            t0 = time.time()
            rec = gczrl.held_out_records(ctx, X, Y, meta, name, radius, calib, verts_of, gcz_of, gt_calib, **kw)
            for lig in LIGAMENTS:
                hctx = ctx  # blend/projection context (zones) for record_error
                for a in ALPHAS:
                    if vname == "rows_fixed_ctx":
                        errs = []
                        for (key, gp, corr, gt, verts) in rec[lig]:
                            e = gczrl.mae(gczrl.blend(ctx_wo[key], gp, corr, a, verts, lig), gt)
                            if e is not None:
                                errs.append(e)
                        err = float(pd.Series(errs).mean()) if errs else float("nan")
                    else:
                        err = gczrl.record_error(hctx, rec[lig], a, lig)
                    grid.append({"variant": vname, "family": fam, "regressor": name, "ligament": lig,
                                 "alpha": a, "mae_mm": err})
            print(f"{vname} {name} {time.time() - t0:.0f}s", flush=True)

grid = pd.DataFrame(grid)
grid.to_csv(OUT / "table3_variants_full_grid.csv", index=False)
summary = {}
for vname in variants:
    g = grid[grid.variant == vname]
    best = g.groupby(["ligament", "family"]).mae_mm.min().unstack()[list(FAMILIES)].reindex(LIGAMENTS)
    best.loc["Mean"] = best.mean()
    best["lowest"] = best.drop(columns=[], errors="ignore").idxmin(axis=1)
    best.round(3).to_csv(OUT / f"table3_{vname}.csv")
    summary[vname] = best.round(3).to_dict()
    print(vname)
    print(best.round(3).to_string())
(OUT / "table3_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
