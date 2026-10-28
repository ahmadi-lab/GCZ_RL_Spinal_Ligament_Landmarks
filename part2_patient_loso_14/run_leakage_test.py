"""Data-leakage test for the patient-level LOSO (Table P2-5).

The S02 fold is recomputed from scratch twice in a temporary folder: once with the true S02
annotations and once with corrupted copies (Gaussian noise, SD 25 mm, on every landmark
coordinate).  If no annotation of the held-out subject reaches training, model selection or
prediction, the GCZ and GCZ-RL predictions for S02 of the two runs must be identical.
Exit code 1 if they are not.

Run after run_part2.py.   Output: outputs/table_P2_5_leakage_test_recomputed.json (the archived
result of the original computing environment, outputs/table_P2_5_leakage_test.json, is kept)
"""
import json
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "code" / "ligament"))
sys.path.insert(0, str(HERE))
import cases  # noqa: E402
import run_part2 as P  # noqa: E402
from data_io import case_code  # noqa: E402
from metrics import to_rows  # noqa: E402

HELD = "Verse561"   # S02
rng = np.random.default_rng(0)
original = cases.landmarks


def corrupted(pat, lv):
    g = original(pat, lv)
    if pat != HELD or g is None:
        return g
    return {k: (None if v is None else v + rng.normal(0, 25, v.shape)) for k, v in g.items()}


def fold_predictions(label_fn, tmp):
    """Recompute the S02 fold from scratch (GCZ geometry and GCZ-RL) with the given label source."""
    cases.landmarks = label_fn
    for sub in ("", "models", "ckpt"):
        (tmp / sub).mkdir(parents=True, exist_ok=True)
    P.OUT, P.CKPT = tmp, tmp / "ckpt"
    gcz_p, rl_p, _, _ = P.run_fold(HELD, fresh=True)
    cases.landmarks = original
    return pd.concat([to_rows(gcz_p, "GCZ", case_code), to_rows(rl_p, "GCZ-RL", case_code)])


def main():
    # Both runs are recomputed from scratch, so the comparison isolates the held-out labels
    # (a checkpointed geometry would otherwise be compared with a recomputed one).
    tmp = HERE / "leakage_test_tmp"
    clean = fold_predictions(original, tmp / "clean")
    new = fold_predictions(corrupted, tmp / "corrupted")
    k = ["method", "case", "ligament", "point_index"]
    m = new.merge(clean, on=k, suffixes=("_new", "_ref"))
    dev = float(np.abs(m[["x_new", "y_new", "z_new"]].to_numpy() - m[["x_ref", "y_ref", "z_ref"]].to_numpy()).max())
    res = {"held_out_subject": "S02", "corruption": "Gaussian noise, SD 25 mm, added to every S02 landmark",
           "recomputed": "GCZ geometry and GCZ-RL training of the S02 fold, from scratch, with clean and with corrupted S02 labels",
           "points_compared": int(len(m)), "points_clean_run": int(len(clean)), "points_corrupted_run": int(len(new)),
           "max_abs_coordinate_difference_mm": dev, "tolerance_mm": 1e-9,
           "predictions_unchanged": bool(dev <= 1e-9 and len(m) == len(clean) == len(new))}
    # table_P2_5_leakage_test.json is the archived result of the original computing environment (424
    # points); this run, with the ICP recomputed here, is written beside it and never overwrites it.
    out = HERE / "outputs"
    (out / "table_P2_5_leakage_test_recomputed.json").write_text(json.dumps(res, indent=2), encoding="utf-8")
    shutil.rmtree(tmp, ignore_errors=True)
    print(json.dumps(res, indent=2))
    if not res["predictions_unchanged"]:
        print("FAIL: S02 predictions depend on S02 annotations")
        sys.exit(1)
    print("PASS: S02 predictions do not depend on S02 annotations")


if __name__ == "__main__":
    main()
