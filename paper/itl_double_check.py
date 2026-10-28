"""Independent double-check of the ITL result (Sections 4.3-4.4, Fig. 7, Supplementary S19).

GCZ-RL lowers the ITL error from 6.262 to 1.796 mm over the cohort and from 8.095 to 0.678 mm on
S03 L4. Because the reduction is large, it is checked here without the evaluation code of the
pipeline:

1. Errors recomputed from the raw annotation files (3D Slicer .mrk.json, read with json) and the
   saved predictions, with an exhaustive one-to-one matching, for the 14 held-out vertebrae and the
   three methods; compared with the metrics file and the reported class means.
2. Position on the mesh of S03 L4 (read with PyVista): each predicted point is a mesh vertex; distance
   of the expert point, the GCZ point, and the GCZ-RL point to the tip of the transverse process
   (the most lateral vertex on its side).
3. Leakage test on the S03 fold: the fold is recomputed from scratch with the true S03 annotations
   and with every S03 annotation replaced by noise (SD 25 mm); the S03 predictions must be identical,
   and the recomputation with the true annotations is compared with the saved predictions.
4. Why the GCZ points miss the tips: the GCZ ITL step (rigid ICP of the atlas, transfer of the atlas
   expert points, snap into the selection region) is recomputed for the 14 vertebrae, and the GCZ error
   is related to how much farther the target's transverse process reaches than the atlas's.
5. The S03 fold without the side constraint.

Usage: python paper/itl_double_check.py [--no-leakage]
Writes paper/outputs/revision/itl_double_check.{json,csv} and itl_gcz_origin.csv; exit code 1 if any
check fails.
"""
from __future__ import annotations

import itertools
import json
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pyvista as pv

REPO = Path(__file__).resolve().parents[1]
P2 = REPO / "part2_patient_loso_14"
OUT = REPO / "paper" / "outputs" / "revision"
sys.path.insert(0, str(REPO / "code" / "ligament"))
from data_io import RAW_ROOT  # noqa: E402  (path of the raw data only)

PATIENT = {"S01": "Verse005", "S02": "Verse561", "S03": "Verse581"}
TOL = 1e-9


def expert_itl(pat, lv):
    """ITL expert points per side, read directly from the .mrk.json files (LPS, mm)."""
    for sub in (f"{pat}_STL", pat):
        for suffix in ("_ligaments_GT", "_ligaments"):
            d = RAW_ROOT / sub / f"{pat}_{lv}{suffix}"
            if d.exists():
                out = {}
                for side in ("L", "R"):
                    files = sorted(d.glob(f"{pat}_{lv}_ITL_{side}*.mrk.json"))
                    pts = [cp["position"] for f in files
                           for cp in json.loads(f.read_text(encoding="utf-8"))["markups"][0]["controlPoints"]
                           if cp.get("positionStatus", "defined") == "defined"]
                    out[side] = np.asarray(pts, float)
                return out
    raise FileNotFoundError(f"no annotations for {pat} {lv}")


def matched(pred, ref):
    """Exhaustive one-to-one matching (min(|P|,|R|) pairs of least total distance); returns the distances."""
    d = np.linalg.norm(pred[:, None] - ref[None], axis=2)
    if len(pred) <= len(ref):
        best = min(itertools.permutations(range(len(ref)), len(pred)), key=lambda c: d[range(len(pred)), c].sum())
        return d[range(len(pred)), best], best
    best = min(itertools.permutations(range(len(pred)), len(ref)), key=lambda c: d[c, range(len(ref))].sum())
    return d[best, range(len(ref))], best


def recompute():
    pred = {"GCZ": pd.read_csv(P2 / "outputs/predictions_gcz.csv"), "GCZ-RL": pd.read_csv(P2 / "outputs/predictions_gczrl.csv"),
            "SLD": pd.read_csv(P2 / "outputs/predictions_sld.csv")}
    abl = pd.read_csv(REPO / "paper/outputs/revision/ablation_per_case.csv")
    pred_proj = abl[(abl.method == "GCZ + projection only (no residual)") & (abl.ligament == "ITL")].set_index("case").mae_mm
    met = pd.read_csv(P2 / "outputs/per_vertebra_per_ligament_metrics.csv")
    met = met[met.ligament == "ITL"].set_index(["method", "case"]).mae_mm
    rows = []
    for case in sorted(met.index.get_level_values("case").unique()):
        pat, lv = case.split("_")
        ref = expert_itl(pat, lv)
        refs = np.vstack([ref["L"], ref["R"]])
        sides = ["L"] * len(ref["L"]) + ["R"] * len(ref["R"])
        r = {"case": case, "subject": [s for s, p in PATIENT.items() if p == pat][0], "level": lv, "n_ref": len(refs)}
        for m, df in pred.items():
            p = df[(df.patient == pat) & (df.level == lv) & (df.ligament == "ITL")].sort_values("point_index")
            dist, _ = matched(p[["x", "y", "z"]].to_numpy(float), refs)
            r[f"{m}_mae_recomputed"] = float(dist.mean())
            r[f"{m}_mae_metrics_file"] = float(met.loc[(m, case)])
            r[f"{m}_n_pred"] = len(p)
        r["projection_only_mae"] = float(pred_proj.loc[case])
        if case == "Verse581_L4":
            r["_refs"], r["_sides"] = refs, sides
        rows.append(r)
    return rows


def mesh_check(rows):
    """S03 L4: predicted points on the mesh and their distance to the transverse-process tips."""
    surf = pv.read(str(RAW_ROOT / "Verse581" / "Verse581_L4.vtk")).extract_surface().compute_normals(
        point_normals=True, cell_normals=False, auto_orient_normals=True)
    v, normals = np.asarray(surf.points, float), np.asarray(surf.point_data["Normals"], float)
    mid = 0.5 * (v[:, 0].min() + v[:, 0].max())
    tip = {"R": v[v[:, 0].argmin()], "L": v[v[:, 0].argmax()]}            # LPS: +x is the patient's left
    r = next(x for x in rows if x["case"] == "Verse581_L4")
    refs, sides = r["_refs"], r["_sides"]
    out = []
    for m, f in (("GCZ", "predictions_gcz.csv"), ("GCZ-RL", "predictions_gczrl.csv")):
        p = pd.read_csv(P2 / "outputs" / f)
        p = p[(p.case == "S03_L4") & (p.ligament == "ITL")].sort_values("point_index")[["x", "y", "z"]].to_numpy(float)
        dist, match = matched(p, refs)
        for i, q in enumerate(p):
            side = "L" if q[0] > mid else "R"
            ref = refs[match[i]]
            nb = normals[np.linalg.norm(v - q, axis=1) < 1.5].mean(axis=0)     # outward normal of the surface around the point
            nb /= np.linalg.norm(nb)
            out.append({"method": m, "point": i, "side_of_point": side, "side_of_matched_reference": sides[match[i]],
                        "surface_faces": "posterior" if nb[1] > 0.8 else ("lateral" if abs(nb[0]) > 0.8 else "other"),
                        "normal_x_left": float(nb[0]), "normal_y_posterior": float(nb[1]), "normal_z_superior": float(nb[2]),
                        "error_mm": float(dist[i]), "distance_to_mesh_vertex_mm": float(np.linalg.norm(v - q, axis=1).min()),
                        "point_to_tip_mm": float(np.linalg.norm(q - tip[side])),
                        "reference_to_tip_mm": float(np.linalg.norm(ref - tip[side])),
                        "point_posterior_of_reference_mm": float(q[1] - ref[1])})
    return out


def leakage_s03():
    sys.path.insert(0, str(P2))
    import run_leakage_test as L  # noqa: E402
    L.HELD = "Verse581"
    tmp = P2 / "leakage_test_tmp_S03"
    try:
        clean = L.fold_predictions(L.original, tmp / "clean")
        noisy = L.fold_predictions(L.corrupted, tmp / "corrupted")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    k = ["method", "case", "ligament", "point_index"]
    xyz = ["x", "y", "z"]
    m = noisy.merge(clean, on=k, suffixes=("_n", "_c"))
    saved = pd.concat([pd.read_csv(P2 / "outputs/predictions_gcz.csv"), pd.read_csv(P2 / "outputs/predictions_gczrl.csv")])
    s = clean.merge(saved, on=k, suffixes=("_c", "_s"))
    itl = s[(s.ligament == "ITL") & (s.case == "S03_L4")]
    return {"points_compared_noisy_vs_true": int(len(m)), "points_true_run": int(len(clean)), "points_noisy_run": int(len(noisy)),
            "max_difference_noisy_vs_true_mm": float(np.abs(m[[c + "_n" for c in xyz]].to_numpy() - m[[c + "_c" for c in xyz]].to_numpy()).max()),
            "points_compared_with_saved": int(len(s)),
            "max_difference_recomputed_vs_saved_mm": float(np.abs(s[[c + "_c" for c in xyz]].to_numpy() - s[[c + "_s" for c in xyz]].to_numpy()).max()),
            "max_difference_recomputed_vs_saved_S03_L4_ITL_mm": float(np.abs(itl[[c + "_c" for c in xyz]].to_numpy() - itl[[c + "_s" for c in xyz]].to_numpy()).max())}


def gcz_origin():
    """Where GCZ puts the ITL points: the atlas expert points are moved by the rigid ICP of the fold onto the
    target and snapped to the nearest vertex of the ITL selection region on their side. Recomputed for the
    14 held-out vertebrae and compared with the saved GCZ points; the lateral reach of the transverse process
    is the distance of the expert ITL point from the vertex centroid in x (LPS)."""
    import pickle
    sys.path.insert(0, str(P2))
    import cases  # noqa: E402
    import gcz  # noqa: E402
    import run_part2 as R  # noqa: E402
    rows = []
    for subj, pat in PATIENT.items():
        geo = pickle.loads((P2 / f"checkpoints/fold_held_{subj}_gcz.pkl").read_bytes())
        train_keys, test_keys = R.keys_of([p for p in PATIENT.values() if p != pat]), R.keys_of([pat])
        thr = gcz.build_thresholds(cases.train_dict(train_keys, "vtk"))
        for i, key in enumerate(test_keys):     # same ICP settings and seeds as run_part2.gcz_predict (seed0=1000)
            v = cases.mesh(*key, "vtk")
            akey = next(k for k in train_keys if R.case_code(*k) == geo["test_atlas"][key].split(" ")[0])
            av, agt, tgt = cases.mesh(*akey, "vtk"), cases.landmarks(*akey), cases.landmarks(*key)
            r, t = gcz.icp(av, v, n_sub=5000, max_iter=150, tol=1e-6, seed=200 + 1000 + i)
            sel = gcz.candidate_vertices(v, "ITL", thr, third_fallback=True)
            proj = gcz.candidate_vertices(v, "ITL", thr, third_fallback=False)
            cx, acx = v[:, 0].mean(), av[:, 0].mean()
            for side in ("R", "L"):
                on = lambda p, c: (p[:, 0] >= c) == (side == "L")  # noqa: E731
                a, ref = agt["ITL"][on(agt["ITL"], acx)][0], tgt["ITL"][on(tgt["ITL"], cx)][0]
                warped = gcz.apply(r, t, a[None])[0]
                half = sel[on(sel, cx)]
                snapped = half[np.linalg.norm(half - warped, axis=1).argmin()]
                saved = geo["test_gcz"][key]["ITL"]
                tip = v[v[:, 0].argmax()] if side == "L" else v[v[:, 0].argmin()]
                ph = proj[on(proj, cx)]
                rows.append({"case": R.case_code(*key), "atlas": R.case_code(*akey), "side": side,
                             "atlas_reach_mm": abs(a[0] - acx), "target_reach_mm": abs(ref[0] - cx),
                             "transferred_point_medial_of_reference_mm": (ref[0] - warped[0]) * (1 if side == "L" else -1),
                             "gcz_error_mm": float(np.linalg.norm(snapped - ref)),
                             "resnapped_equals_saved_gcz": bool(np.abs(saved - snapped).sum(axis=1).min() < 1e-6),
                             "selection_region_vertices": len(sel), "projection_region_vertices_side": len(ph),
                             "projection_region_max_to_tip_mm": float(np.linalg.norm(ph - tip, axis=1).max()) if len(ph) else np.nan,
                             "reference_to_tip_mm": float(np.linalg.norm(ref - tip))})
    d = pd.DataFrame(rows)
    d["reach_difference_mm"] = d.target_reach_mm - d.atlas_reach_mm      # > 0: the atlas process is shorter
    d.to_csv(OUT / "itl_gcz_origin.csv", index=False)
    wide = d[d.selection_region_vertices > 1000]
    return d, {"points": len(d), "resnapped_equal_saved": int(d.resnapped_equals_saved_gcz.sum()), "wide_region_points": len(wide),
               "r_error_vs_shortfall_wide": float(np.corrcoef(wide.gcz_error_mm, wide.reach_difference_mm.clip(lower=0))[0, 1]),
               "r_error_vs_abs_difference_wide": float(np.corrcoef(wide.gcz_error_mm, wide.reach_difference_mm.abs())[0, 1]),
               "S03_L4": d[d.case == "S03_L4"].set_index("side").round(4).to_dict(orient="index"),
               "S01_L1": d[d.case == "S01_L1"].set_index("side").round(4).to_dict(orient="index"),
               "S03_L1": d[d.case == "S03_L1"].set_index("side").round(4).to_dict(orient="index")}


def without_side_constraint_s03(rows):
    """The S03 fold re-run with the side constraint switched off (stored geometry, nested selection included)."""
    import tempfile
    sys.path.insert(0, str(P2))
    import gczrl  # noqa: E402
    import run_part2  # noqa: E402
    keep_lock, keep_out = set(gczrl.SIDE_LOCKED), run_part2.OUT
    with tempfile.TemporaryDirectory() as tmp:
        run_part2.OUT = Path(tmp)
        (Path(tmp) / "models").mkdir()
        gczrl.SIDE_LOCKED = set()
        try:
            _, pred, _, _ = run_part2.run_fold("Verse581", fresh=False)
        finally:
            gczrl.SIDE_LOCKED, run_part2.OUT = keep_lock, keep_out
    out = {}
    for r in rows:
        pat, lv = r["case"].split("_")
        if pat == "Verse581":
            ref = expert_itl(pat, lv)
            out[r["case"]] = {"with_constraint": round(r["GCZ-RL_mae_recomputed"], 3),
                              "without_constraint": round(float(matched(pred[(pat, lv)]["ITL"], np.vstack([ref["L"], ref["R"]]))[0].mean()), 3)}
    return out


def main():
    rows = recompute()
    tab = pd.DataFrame([{k: v for k, v in r.items() if not k.startswith("_")} for r in rows])
    tab.to_csv(OUT / "itl_double_check.csv", index=False)
    pts = mesh_check(rows)
    nums = json.loads((REPO / "paper/outputs/paper_numbers.json").read_text())
    res = {"per_vertebra": "itl_double_check.csv", "S03_L4_points": pts}
    checks = {}
    for m in ("GCZ", "GCZ-RL", "SLD"):
        dev = float((tab[f"{m}_mae_recomputed"] - tab[f"{m}_mae_metrics_file"]).abs().max())
        res[f"{m}_ITL_class_mean_recomputed"] = float(tab[f"{m}_mae_recomputed"].mean())
        checks[f"{m}: recomputed from raw files = metrics file (max diff {dev:.1e} mm)"] = dev < TOL
    s = tab.set_index("case").loc["Verse581_L4"]
    res["S03_L4"] = {m: float(s[f"{m}_mae_recomputed"]) for m in ("GCZ", "GCZ-RL", "SLD")} | {"projection_only": float(s.projection_only_mae)}
    checks["class means as reported (6.262 / 1.796 / 0.873 mm)"] = (
        [round(res[f"{m}_ITL_class_mean_recomputed"], 3) for m in ("GCZ", "GCZ-RL", "SLD")] == [6.262, 1.796, 0.873])
    checks["S03 L4 as reported (8.095 / 0.678 mm)"] = (round(res["S03_L4"]["GCZ"], 3), round(res["S03_L4"]["GCZ-RL"], 3)) == (8.095, 0.678)
    checks["every S03 L4 point is a mesh vertex"] = all(p["distance_to_mesh_vertex_mm"] < 1e-3 for p in pts)
    checks["GCZ-RL points matched to the reference of their own side"] = all(
        p["side_of_point"] == p["side_of_matched_reference"] for p in pts if p["method"] == "GCZ-RL")
    checks["GCZ points on the posterior surface of the processes, GCZ-RL points at the lateral tips (Fig. 7 caption)"] = all(
        p["surface_faces"] == ("posterior" if p["method"] == "GCZ" else "lateral") for p in pts)
    red = (tab["GCZ_mae_recomputed"] - tab["GCZ-RL_mae_recomputed"])
    res["ITL_reduction_per_vertebra_mm"] = {c: round(float(x), 3) for c, x in zip(tab.case, red)}
    res["S03_L4_rank_of_ITL_reduction"] = int((red > red[tab.case == "Verse581_L4"].iloc[0]).sum() + 1)
    res["vertebrae_lower_with_GCZRL"] = int((red > 0).sum())
    res["GCZ_ITL_error_range_mm"] = [float(tab["GCZ_mae_recomputed"].min()), float(tab["GCZ_mae_recomputed"].max())]
    _, origin = gcz_origin()
    res["gcz_origin"] = origin
    checks["GCZ ITL step recomputed (ICP, transfer, snap) = saved GCZ points"] = origin["resnapped_equal_saved"] == origin["points"]
    side = without_side_constraint_s03(rows)
    res["S03_ITL_without_side_constraint"] = side
    checks["S03 L4 unchanged without the side constraint"] = side["Verse581_L4"]["without_constraint"] == side["Verse581_L4"]["with_constraint"]
    if "--no-leakage" not in sys.argv:
        leak = leakage_s03()
        res["leakage_test_S03"] = leak
        checks["S03 predictions identical with noisy S03 annotations"] = (
            leak["max_difference_noisy_vs_true_mm"] <= TOL and leak["points_compared_noisy_vs_true"] == leak["points_true_run"] == leak["points_noisy_run"])
        checks["S03 fold recomputed from scratch = saved predictions"] = leak["max_difference_recomputed_vs_saved_mm"] <= TOL
    res["checks"] = checks
    (OUT / "itl_double_check.json").write_text(json.dumps(res, indent=2), encoding="utf-8")
    print(tab[["case", "GCZ_mae_recomputed", "projection_only_mae", "GCZ-RL_mae_recomputed", "SLD_mae_recomputed"]].round(3).to_string(index=False))
    print(pd.DataFrame(pts).round(3).to_string(index=False))
    for c, ok in checks.items():
        print(("PASS  " if ok else "FAIL  ") + c)
    sys.exit(0 if all(checks.values()) else 1)


if __name__ == "__main__":
    main()
