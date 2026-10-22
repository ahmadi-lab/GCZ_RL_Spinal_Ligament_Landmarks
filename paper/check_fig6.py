"""Check Fig. 6 (held-out vertebra S03 L4, primary LOSO) against the saved outputs, and
produce the numbers quoted in its caption.

1. the points drawn by fig6_qualitative_loso.py are exactly the saved predictions of each method;
2. every drawn point lies on the vertebral surface (also confirms the SLD frame conversion);
3. errors recomputed from those points equal per_vertebra_per_ligament_metrics.csv;
4. the caption numbers: per-method error on this vertebra, its rank among the 14 held-out
   vertebrae by GCZ -> GCZ-RL improvement.

Writes paper/outputs/fig6_check_S03_L4.csv and paper/outputs/fig6_vertebra_ranking.csv;
exit code 1 if any check fails.   Usage:  python paper/check_fig6.py
"""
import sys
from pathlib import Path

import pandas as pd
from scipy.spatial import cKDTree

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "code" / "ligament"))
import cases  # noqa: E402
from data_io import LIGAMENTS  # noqa: E402
from metrics import pair_metrics  # noqa: E402

P2 = REPO / "part2_patient_loso_14" / "outputs"
OUT = REPO / "paper" / "outputs"
PAT, LV, CASE = "Verse581", "L4", "S03_L4"
METHODS = [("GCZ", "predictions_gcz.csv"), ("GCZ-RL", "predictions_gczrl.csv"), ("SLD", "predictions_sld.csv")]
SURFACE_TOL_MM = 1.0   # SLD landmarks are snapped to the surface, not necessarily to a vertex


def main():
    v = cases.mesh(PAT, LV, "vtk")
    kd = cKDTree(v)
    gt = cases.landmarks(PAT, LV)
    saved = pd.read_csv(P2 / "per_vertebra_per_ligament_metrics.csv")
    saved_case = saved[saved.case == f"{PAT}_{LV}"].set_index(["method", "ligament"])

    rows = []
    for m, fn in METHODS:
        d = pd.read_csv(P2 / fn)
        d = d[d.case == CASE]
        for lig in LIGAMENTS:
            p = d[d.ligament == lig].sort_values("point_index")[["x", "y", "z"]].to_numpy(float)
            pm = pair_metrics(p, gt.get(lig))
            rows.append({"method": m, "ligament": lig, "n_pred": len(p), "n_ref": pm["n_ref"],
                         "mae_recomputed": round(pm["mae_mm"], 3),
                         "mae_saved": round(float(saved_case.loc[(m, lig), "mae_mm"]), 3),
                         "max_dist_to_surface_mm": round(float(kd.query(p)[0].max()), 4)})
    t = pd.DataFrame(rows)
    OUT.mkdir(parents=True, exist_ok=True)
    t.to_csv(OUT / "fig6_check_S03_L4.csv", index=False)
    print(t.to_string(index=False))

    vm = saved.pivot_table(index="case", columns="method", values="mae_mm", aggfunc="mean")
    vm["improvement_GCZ_to_GCZRL"] = vm["GCZ"] - vm["GCZ-RL"]
    vm = vm.sort_values("improvement_GCZ_to_GCZRL", ascending=False).round(3)
    vm.to_csv(OUT / "fig6_vertebra_ranking.csv")
    rank = list(vm.index).index(f"{PAT}_{LV}") + 1
    macro = t.groupby("method").mae_recomputed.mean().round(3).to_dict()

    ok_mae = bool((t.mae_recomputed - t.mae_saved).abs().max() < 1e-3)
    ok_surface = bool(t.max_dist_to_surface_mm.max() <= SURFACE_TOL_MM)
    print(f"\nS03 L4 macro MAE: GCZ {macro['GCZ']}, GCZ-RL {macro['GCZ-RL']}, SLD {macro['SLD']} mm")
    print(f"rank of S03 L4 by GCZ -> GCZ-RL improvement: {rank} of {len(vm)}")
    print(f"{'PASS' if ok_mae else 'FAIL'}  recomputed errors equal the saved metrics")
    print(f"{'PASS' if ok_surface else 'FAIL'}  every drawn point lies on the surface "
          f"(max {t.max_dist_to_surface_mm.max()} mm)")
    expected = {"GCZ": 3.263, "GCZ-RL": 1.972, "SLD": 2.771}
    ok_caption = rank == 5 and all(abs(macro[k] - val) < 1e-3 for k, val in expected.items())
    print(f"{'PASS' if ok_caption else 'FAIL'}  caption numbers (3.263 / 1.972 / 2.771 mm, rank 5 of 14)")
    sys.exit(0 if ok_mae and ok_surface and ok_caption else 1)


if __name__ == "__main__":
    main()
