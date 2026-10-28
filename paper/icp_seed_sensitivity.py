"""Sensitivity of the primary evaluation to the random seeds of the ICP atlas selection.

The complete nested LOSO procedure of part2_patient_loso_14/run_part2.py is re-run with the ICP
geometry recomputed from scratch (atlas selection, alignment, training-side GCZ predictions), with
the random seeds of the ICP point sampling shifted by a fixed offset. Nothing else changes: same
code, same settings, same inner model selection, no use of held-out annotations before the final
scoring. Offset 0 uses the seeds of the method; the stored checkpoints are neither read nor
written (a temporary checkpoint folder is used). SLD does not use ICP and is not re-run.

Usage:  python paper/icp_seed_sensitivity.py run 0 1 2 ...   (one JSON per offset)
        python paper/icp_seed_sensitivity.py summary
Outputs: paper/outputs/icp_seed_sensitivity/offset_<k>.json, summary.csv, summary.json, tab_icp_seeds.tex
(the figure, Fig. 8, is drawn by paper/figures_main.py)
"""
from __future__ import annotations

import json
import sys
import tempfile
import time
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
P2 = REPO / "part2_patient_loso_14"
OUT = REPO / "paper" / "outputs" / "icp_seed_sensitivity"
sys.path.insert(0, str(REPO / "code" / "ligament"))
sys.path.insert(0, str(P2))
import cases  # noqa: E402
import run_part2  # noqa: E402
from data_io import LIGAMENTS, SUBJECT_CODE  # noqa: E402
from metrics import evaluate  # noqa: E402

SEED_STEP = 10_000          # offset k shifts every ICP seed stream of gcz_predict by k * SEED_STEP
SLD_MACRO_MAE = 2.224       # SLD has no ICP step; value of the primary evaluation


def scores(pred, gt):
    df = evaluate(pred, gt, "m", subject_of=lambda k: SUBJECT_CODE[k[0]])
    ok = df[~df.missing_prediction]
    macro = lambda d: float(d.groupby("ligament").mae_mm.mean().reindex(LIGAMENTS).mean())  # noqa: E731
    return {"macro_mae_mm": macro(ok),
            "pass5_case_avg_pct": float(ok.groupby("ligament").pass5.mean().reindex(LIGAMENTS).mean()),
            "ref_to_pred_mm": float(ok.ref_to_pred_mean_mm.mean()),
            "by_subject": {s: macro(x) for s, x in ok.groupby("subject")},
            "by_class": {k: float(v) for k, v in ok.groupby("ligament").mae_mm.mean().items()}}


def run(offset: int):
    original = run_part2.gcz_predict

    def shifted(targets, atlas_keys, zone_keys, seed0):
        return original(targets, atlas_keys, zone_keys, seed0 + offset * SEED_STEP)

    stored = {}
    for f in json.loads((P2 / "outputs" / "fold_details.json").read_text()):
        stored.update(f["test_atlas_choice"])
    t0 = time.time()
    with tempfile.TemporaryDirectory() as tmp:
        run_part2.OUT, run_part2.CKPT = Path(tmp), Path(tmp) / "ckpt"
        for d in (run_part2.OUT / "models", run_part2.CKPT):
            d.mkdir(parents=True, exist_ok=True)
        run_part2.gcz_predict = shifted
        pred_gcz, pred_rl, atlas, selection = {}, {}, {}, {}
        for subj in run_part2.COHORT:
            g, r, info, _ = run_part2.run_fold(subj, fresh=True)
            pred_gcz.update(g)
            pred_rl.update(r)
            atlas.update(info["test_atlas_choice"])
            selection[info["held_out_subject"]] = {l: [v["extra_trees"], v["alpha"]] for l, v in info["selection"].items()}
        run_part2.gcz_predict = original
    gt = {k: cases.landmarks(*k) for k in pred_rl}
    res = {"offset": offset, "seconds": round(time.time() - t0),
           "GCZ": scores(pred_gcz, gt), "GCZ-RL": scores(pred_rl, gt),
           "atlas": atlas, "atlas_changed_vs_stored": sorted(k for k in atlas if atlas[k].split()[0] != stored[k].split()[0]),
           "selection": selection}
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / f"offset_{offset}.json").write_text(json.dumps(res, indent=2), encoding="utf-8")
    print(f"offset {offset}: GCZ {res['GCZ']['macro_mae_mm']:.3f}  GCZ-RL {res['GCZ-RL']['macro_mae_mm']:.3f}  "
          f"atlas changed {len(res['atlas_changed_vs_stored'])}/14  ({res['seconds']} s)", flush=True)


def summary():
    rows = []
    for p in sorted(OUT.glob("offset_*.json"), key=lambda p: int(p.stem.split("_")[1])):
        d = json.loads(p.read_text())
        g, r = d["GCZ"], d["GCZ-RL"]
        rows.append({"offset": d["offset"], "GCZ_mm": g["macro_mae_mm"], "GCZRL_mm": r["macro_mae_mm"],
                     "reduction_pct": 100 * (g["macro_mae_mm"] - r["macro_mae_mm"]) / g["macro_mae_mm"],
                     "GCZ_pass5": g["pass5_case_avg_pct"], "GCZRL_pass5": r["pass5_case_avg_pct"],
                     "GCZRL_R2P_mm": r["ref_to_pred_mm"],
                     "subjects_improved": sum(r["by_subject"][s] < g["by_subject"][s] for s in g["by_subject"]),
                     "GCZRL_below_SLD": r["macro_mae_mm"] < SLD_MACRO_MAE,
                     "atlas_changed": len(d["atlas_changed_vs_stored"]),
                     **{f"GCZRL_{s}": v for s, v in r["by_subject"].items()}})
    df = pd.DataFrame(rows)
    df.round(3).to_csv(OUT / "summary.csv", index=False)
    stats = {c: {"mean": float(df[c].mean()), "sd": float(df[c].std(ddof=1)) if len(df) > 1 else 0.0,
                 "min": float(df[c].min()), "max": float(df[c].max()), "median": float(df[c].median())}
             for c in ("GCZ_mm", "GCZRL_mm", "reduction_pct", "GCZRL_pass5", "GCZRL_R2P_mm", "atlas_changed")}
    stats["runs"] = int(len(df))
    stats["runs_all_subjects_improved"] = int((df.subjects_improved == 3).sum())
    stats["runs_GCZRL_below_SLD"] = int(df.GCZRL_below_SLD.sum())
    (OUT / "summary.json").write_text(json.dumps(stats, indent=2), encoding="utf-8")
    lines = [f"{int(r.offset)} & {r.GCZ_mm:.3f} & {r.GCZRL_mm:.3f} & {r.reduction_pct:.1f} & {r.GCZRL_pass5:.1f} & "
             f"{r.GCZRL_S01:.3f} & {r.GCZRL_S02:.3f} & {r.GCZRL_S03:.3f} & {int(r.subjects_improved)}/3 & "
             f"{int(r.atlas_changed)}/14 \\\\" for r in df.itertuples()]
    lines += ["\\midrule",
              f"Mean $\\pm$ SD & {df.GCZ_mm.mean():.3f} $\\pm$ {df.GCZ_mm.std(ddof=1):.3f} & "
              f"{df.GCZRL_mm.mean():.3f} $\\pm$ {df.GCZRL_mm.std(ddof=1):.3f} & {df.reduction_pct.mean():.1f} & "
              f"{df.GCZRL_pass5.mean():.1f} & {df.GCZRL_S01.mean():.3f} & {df.GCZRL_S02.mean():.3f} & "
              f"{df.GCZRL_S03.mean():.3f} & {int((df.subjects_improved == 3).sum())}/{len(df)} runs & "
              f"{df.atlas_changed.mean():.1f} \\\\"]
    (OUT / "tab_icp_seeds.tex").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(df.round(3).to_string(index=False))
    print(json.dumps(stats, indent=1))


if __name__ == "__main__":
    if sys.argv[1] == "run":
        for k in sys.argv[2:]:
            run(int(k))
    else:
        summary()
