"""Reported-versus-reproduced table for the main numerical claims of the manuscript.

Every reported value below is copied from the manuscript text or tables. The reproduced value is
read from the outputs that the run scripts regenerate, so this script must run after them:

  level A  metrics recomputed from the saved predictions and the annotations
           (paper/make_paper_numbers.py, code/verify_all.py)
  level B  end-to-end re-run with the stored ICP alignments and SLD output; all Extra Trees models
           are retrained with their fixed seeds (part2 run_part2.py, part1 run_part1.py,
           paper/table3_regressor_comparison.py)
  level C  ICP geometry recomputed from scratch (leakage test of the S02 fold; archived
           second-environment run in paper/outputs/sensitivity_fresh_icp, not re-executed here)

Writes paper/outputs/revision/reproducibility_table.csv and .tex; exit code 1 on any mismatch.
Usage: python paper/reproducibility_check.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[1]
P1, P2, PAPER = REPO / "part1_dissertation_57" / "outputs", REPO / "part2_patient_loso_14" / "outputs", REPO / "paper" / "outputs"
OUT = PAPER / "revision"
TOL = 0.0005


def main():
    nums = json.loads((PAPER / "paper_numbers.json").read_text())
    overall = pd.read_csv(P2 / "table_P2_2_overall.csv").set_index("method")
    lig = pd.read_csv(P2 / "table_P2_3_per_ligament.csv").set_index(["method", "ligament"]).MAE_mm
    subj = pd.read_csv(P2 / "table_P2_4_per_subject.csv").set_index(["method", "held_out_subject"]).MAE_mm
    leak = json.loads((P2 / "table_P2_5_leakage_test.json").read_text())
    folds = json.loads((P2 / "fold_details.json").read_text())
    dev = pd.read_csv(P1 / "table_4_4_heldout_results.csv")
    t3 = pd.read_csv(PAPER / "table3_rows_fixed_ctx.csv")
    fresh = pd.read_csv(PAPER / "sensitivity_fresh_icp" / "table_P2_2_overall.csv").set_index("method")
    g, rl, sld = "GCZ", "GCZ-RL", "SLD"
    red = lambda a, b: round(100 * (a - b) / a, 1)  # noqa: E731

    rows = [
        # claim, reported, reproduced, level, artifacts
        ("Primary cohort: vertebrae", 14, sum(len(f["test_vertebrae"]) for f in folds), "B", "fold_details.json"),
        ("Primary cohort: subjects", 3, len(folds), "B", "fold_details.json"),
        ("GCZ macro MAE (mm)", 2.879, overall.loc[g, "MAE_mm"], "B", "table_P2_2_overall.csv"),
        ("GCZ-RL macro MAE (mm)", 2.150, overall.loc[rl, "MAE_mm"], "B", "table_P2_2_overall.csv"),
        ("SLD macro MAE (mm)", 2.224, overall.loc[sld, "MAE_mm"], "B", "table_P2_2_overall.csv"),
        ("GCZ to GCZ-RL error reduction (%)", 25.3, red(overall.loc[g, "MAE_mm"], overall.loc[rl, "MAE_mm"]), "B", "table_P2_2_overall.csv"),
        ("GCZ Pass@5 (%)", 86.0, nums[g]["pass5_case_avg"], "A", "paper_numbers.json"),
        ("GCZ-RL Pass@5 (%)", 92.7, nums[rl]["pass5_case_avg"], "A", "paper_numbers.json"),
        ("SLD Pass@5 (%)", 90.4, nums[sld]["pass5_case_avg"], "A", "paper_numbers.json"),
        ("GCZ reference points within 5 mm (%)", 72.1, nums[g]["ref_within_5mm_pct"], "A", "paper_numbers.json"),
        ("GCZ-RL reference points within 5 mm (%)", 78.4, nums[rl]["ref_within_5mm_pct"], "A", "paper_numbers.json"),
        ("SLD reference points within 5 mm (%)", 77.6, nums[sld]["ref_within_5mm_pct"], "A", "paper_numbers.json"),
        ("GCZ-RL R-to-P closest point (mm)", 3.297, nums[rl]["closest_point_ref_to_pred"], "A", "paper_numbers.json"),
        ("SLD R-to-P closest point (mm)", 4.089, nums[sld]["closest_point_ref_to_pred"], "A", "paper_numbers.json"),
        ("GCZ R-to-P closest point (mm)", 4.024, nums[g]["closest_point_ref_to_pred"], "A", "paper_numbers.json"),
        ("SLD P-to-R closest point (mm)", 1.787, nums[sld]["closest_point_pred_to_ref"], "A", "paper_numbers.json"),
        ("GCZ-RL P-to-R closest point (mm)", 2.140, nums[rl]["closest_point_pred_to_ref"], "A", "paper_numbers.json"),
        ("Reference landmarks inside candidate regions (%)", 90.0, nums["zones_inside_pct"], "A", "paper_numbers.json"),
        ("GCZ-RL minus GCZ (mm)", -0.729, nums["GCZRL_vs_GCZ"]["macro_mae_difference_mm"], "A", "paper_numbers.json"),
        ("  95% bootstrap CI, lower (mm)", -1.078, nums["GCZRL_vs_GCZ"]["bootstrap_ci95_mm"][0], "A", "paper_numbers.json"),
        ("  95% bootstrap CI, upper (mm)", -0.381, nums["GCZRL_vs_GCZ"]["bootstrap_ci95_mm"][1], "A", "paper_numbers.json"),
        ("  Wilcoxon p (reported as < 0.001)", 0.001, nums["GCZRL_vs_GCZ"]["wilcoxon_p"], "A", "paper_numbers.json"),
        ("GCZ-RL minus SLD (mm)", -0.074, nums["GCZRL_vs_SLD"]["macro_mae_difference_mm"], "A", "paper_numbers.json"),
        ("  Wilcoxon p", 0.35, round(nums["GCZRL_vs_SLD"]["wilcoxon_p"], 2), "A", "paper_numbers.json"),
        ("Vertebrae improved by GCZ-RL over GCZ", 13, nums["vertebrae_improved_GCZRL_vs_GCZ"], "A", "paper_numbers.json"),
    ]
    # per-subject reductions from unrounded per-case errors (the submitted draft computed them from rounded means)
    cases = pd.read_csv(P2 / "per_vertebra_per_ligament_metrics.csv")
    cases = cases[~cases.missing_prediction]
    exact = {(m, s): x.groupby("ligament").mae_mm.mean().mean() for (m, s), x in cases.groupby(["method", "subject"])}
    for s, (a, b, pct) in {"S01": (3.784, 2.315, 38.8), "S02": (2.208, 1.986, 10.1), "S03": (2.512, 2.117, 15.7)}.items():
        rows += [(f"{s}: GCZ macro MAE (mm)", a, subj[(g, s)], "B", "table_P2_4_per_subject.csv"),
                 (f"{s}: GCZ-RL macro MAE (mm)", b, subj[(rl, s)], "B", "table_P2_4_per_subject.csv"),
                 (f"{s}: reduction (%)", pct, red(exact[(g, s)], exact[(rl, s)]), "B", "per_vertebra_per_ligament_metrics.csv")]
    table4 = {"ALL": (2.605, 1.658, 2.382), "PLL": (1.863, 2.445, 1.863), "LF": (1.918, 4.997, 1.912),
              "ISL": (2.291, 2.039, 2.226), "SSL": (2.917, 1.968, 2.692), "ITL": (6.262, 0.873, 1.796),
              "CL": (2.302, 1.591, 2.180)}
    for c, vals in table4.items():
        for m, v in zip((g, sld, rl), vals):
            rows.append((f"{c}: {m} MAE (mm)", v, lig[(m, c)], "B", "table_P2_3_per_ligament.csv"))
    rows += [("Leakage test: max coordinate difference (mm)", 0.0, leak["max_abs_coordinate_difference_mm"], "C",
              "table_P2_5_leakage_test.json"),
             ("Leakage test: points compared", 424, leak["points_compared"], "C", "table_P2_5_leakage_test.json")]
    dev = dev.set_index(["method", "ligament"])
    for name, mae, p5 in (("GCZ", 3.708, 79.8), ("GCZ-RL", 2.946, 88.1)):
        rows += [(f"Development: {name} macro MAE (mm)", mae, dev.loc[(name, "Mean"), "MAE_mm"], "B",
                  "table_4_4_heldout_results.csv"),
                 (f"Development: {name} Pass@5 (%)", p5, round(dev.loc[(name, "Mean"), "Pass5_percent"], 1), "B",
                  "table_4_4_heldout_results.csv")]
    t3 = t3.set_index("ligament")
    for fam, rep in (("Extra trees", 2.037), ("Random forest", 2.307), ("Gradient boosting", 2.337)):
        rows.append((f"Calibration comparison: {fam} mean MAE (mm)", rep, t3.loc["Mean", fam], "B",
                     "table3_rows_fixed_ctx.csv"))
    rows += [("Second environment: GCZ macro MAE (mm)", 3.208, fresh.loc[g, "MAE_mm"], "C (archived)",
              "sensitivity_fresh_icp/table_P2_2_overall.csv"),
             ("Second environment: GCZ-RL macro MAE (mm)", 2.418, fresh.loc[rl, "MAE_mm"], "C (archived)",
              "sensitivity_fresh_icp/table_P2_2_overall.csv")]

    seed0 = PAPER / "icp_seed_sensitivity" / "offset_0.json"
    if seed0.exists():  # ICP recomputed here with the method's seeds (executed, not archived)
        s0 = json.loads(seed0.read_text())
        rows += [("Recomputed ICP, seed offset 0: GCZ macro MAE (mm)", 3.208, s0["GCZ"]["macro_mae_mm"], "C",
                  "icp_seed_sensitivity/offset_0.json"),
                 ("Recomputed ICP, seed offset 0: GCZ-RL macro MAE (mm)", 2.418, s0["GCZ-RL"]["macro_mae_mm"], "C",
                  "icp_seed_sensitivity/offset_0.json")]
    df = pd.DataFrame(rows, columns=["claim", "reported", "reproduced", "level", "artifact"])
    df["reproduced"] = df.reproduced.astype(float)
    df["abs_difference"] = (df.reported - df.reproduced).abs().round(4)
    # values are reported to the precision shown; a difference below half a unit of the last digit is a match
    prec = df.reported.map(lambda x: 10 ** -len(repr(float(x)).split(".")[1]) if "." in repr(float(x)) else 1)
    df["status"] = ["verified" if d <= max(p / 2, TOL) + 1e-9 else "MISMATCH" for d, p in zip(df.abs_difference, prec)]
    upper = df.claim.str.contains("reported as <")  # an upper bound is verified when the value lies below it
    df.loc[upper, "abs_difference"] = 0.0
    df.loc[upper, "status"] = ["verified" if v < b else "MISMATCH" for v, b in zip(df.reproduced[upper], df.reported[upper])]
    # rounding corrections made in the revised manuscript (reported value of the submitted draft)
    corrected = {"S02: reduction (%)": "corrected to 10.0 (rounding)"}
    df.loc[df.claim.isin(corrected) & (df.status == "MISMATCH"), "status"] = df.claim.map(corrected)
    OUT.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT / "reproducibility_table.csv", index=False)

    def shown(v, claim):
        if "reported as <" in claim:
            return f"{v:.1e}" if v < 0.001 else f"{v:g}"
        if "(mm)" in claim:
            return f"{v:.3f}"
        if "(%)" in claim:
            return f"{v:.1f}"
        return f"{v:.2f}" if "p" in claim.split() else f"{v:g}"
    with open(OUT / "reproducibility_table.tex", "w", encoding="utf-8") as fh:
        for _, r in df.iterrows():
            claim = r.claim.replace("%", r"\%").replace("&", r"\&").replace("<", "$<$").replace("  ", r"\quad ")
            rep = "$<$0.001" if "reported as <" in r.claim else shown(r.reported, r.claim)
            art = r.artifact.replace("_", "\\_\\allowbreak{}").replace("/", "/\\allowbreak{}")
            fh.write(f"{claim} & {rep} & {shown(r.reproduced, r.claim)} & {r.abs_difference:g} & {r.level} & "
                     f"\\texttt{{{art}}} & {r.status} \\\\\n")
    bad = df[df.status == "MISMATCH"]
    print(df.to_string(index=False))
    print(f"\n{len(df) - len(bad)}/{len(df)} reported values verified "
          f"({int(df.status.str.startswith('corrected').sum())} corrected in the revision)")
    sys.exit(1 if len(bad) else 0)


if __name__ == "__main__":
    main()
