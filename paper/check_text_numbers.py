"""Check the numbers that the revision added to the figure captions and to Sections 4.3-4.4, S13-S19
of the manuscript against the regenerated outputs (the main claims are checked by
paper/reproducibility_check.py).

Each claim is a passage built from the regenerated values with the formatting of the manuscript; it
must occur verbatim in manuscript/main.tex, supplementary.tex, or the generated table rows. A value that changed, or
a number in the text that no output supports, makes the check fail.

Usage: python paper/check_text_numbers.py         Exit code 1 if any passage is missing.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
P2 = REPO / "part2_patient_loso_14" / "outputs"
OUT = REPO / "paper" / "outputs"
TEX = "".join(p.read_text(encoding="utf-8") for p in [REPO / "manuscript" / "main.tex",
                                                      REPO / "manuscript" / "supplementary.tex",
                                                      *sorted((REPO / "manuscript" / "tables").glob("*.tex"))])


def f(x, d=3):
    return f"{x:.{d}f}"


def main():
    nums = json.loads((OUT / "paper_numbers.json").read_text())
    rev = json.loads((OUT / "revision" / "revision_numbers.json").read_text())
    seeds = json.loads((OUT / "icp_seed_sensitivity" / "summary.json").read_text())
    trace = json.loads((OUT / "revision" / "itl_trace_S03_L4.json").read_text())
    dc = json.loads((OUT / "revision" / "itl_double_check.json").read_text())
    dct = pd.read_csv(OUT / "revision" / "itl_double_check.csv").set_index("case")
    tip = {(p["method"], p["side_of_point"]): p for p in dc["S03_L4_points"]}
    og = dc["gcz_origin"]
    words = {2: "two", 3: "three", 4: "four", 5: "five", 6: "six"}
    itl_red = dct["GCZ_mae_recomputed"] - dct["GCZ-RL_mae_recomputed"]
    bigger = dct[itl_red > itl_red["Verse581_L4"]].sort_values(["subject", "level"])
    lv = [int(x[1:]) for x in bigger.level]
    bigger_names = (f"{bigger.subject.iloc[0]} L{lv[0]}--L{lv[-1]}"
                    if bigger.subject.nunique() == 1 and lv == list(range(lv[0], lv[-1] + 1))
                    else ", ".join(f"{s} {l}" for s, l in zip(bigger.subject, bigger.level)))
    tr = pd.read_csv(OUT / "revision" / "itl_trace_S03_L4.csv")
    margin = json.loads((OUT / "margin_verification" / "margin_summary.json").read_text())
    zn = pd.read_csv(OUT / "margin_verification" / "margin_zn_by_vertebra.csv")
    reg = pd.read_csv(OUT / "region_verification" / "region_verification_by_class.csv").set_index("ligament")
    v7 = pd.read_csv(OUT / "figures" / "fig7_landmark_verification_S03_L4.csv")
    rank = pd.read_csv(OUT / "fig6_vertebra_ranking.csv")
    cases = pd.read_csv(P2 / "per_vertebra_per_ligament_metrics.csv")
    cases = cases[~cases.missing_prediction]
    s34 = cases[cases.case == "Verse581_L4"].pivot(index="ligament", columns="method", values="mae_mm")
    vert = {(m, c): x.mae_mm.mean() for (m, c), x in cases.groupby(["method", "case"])}
    subj = {(m, s): x.groupby("ligament").mae_mm.mean().mean() for (m, s), x in cases.groupby(["method", "subject"])}
    red = {s: 100 * (subj[("GCZ", s)] - subj[("GCZ-RL", s)]) / subj[("GCZ", s)] for s in ("S01", "S02", "S03")}
    stage = lambda name, side=None: (tr[(tr.stage == name) & (tr.side == side)].error_mm.iloc[0] if side else  # noqa: E731
                                     trace["class_mae_mm"][name])
    itl = v7[v7.ligament == "ITL"].set_index("side")
    cl = v7[v7.ligament == "CL"]
    cl_s = cl[cl.gczrl_ref_part == "CL_S_R"].sort_values("point_index")
    cl_big = cl[cl.correction_mm > 1.0].set_index("side")
    g, rl = "GCZ", "GCZ-RL"
    s03l1 = (vert[(g, "Verse581_L1")], vert[(rl, "Verse581_L1")])
    pos = rank.reset_index(drop=True)
    rank_s03l4 = int(pos.index[pos.case == "Verse581_L4"][0]) + 1
    ords = {1: "first", 2: "second", 3: "third", 4: "fourth", 5: "fifth", 6: "sixth"}
    proj = rev["GCZRL_vs_projection_only"]
    abl = {r["configuration"]: r["macro_mae_mm"] for r in rev["ablation"]}
    lo_i = zn[(zn.ligament == "ALL") & (zn.part == "I")].z_n_mean
    lo_s = zn[(zn.ligament == "ALL") & (zn.part == "S")].z_n_mean
    claims = [
        # Fig. 6 and Section 4.3
        ("Fig. 6: per-subject values", f"S01 {f(subj[(g, 'S01')])} to {f(subj[(rl, 'S01')])}~mm, S02 {f(subj[(g, 'S02')])} to "
                                       f"{f(subj[(rl, 'S02')])}~mm, S03 {f(subj[(g, 'S03')])} to {f(subj[(rl, 'S03')])}~mm"),
        ("Sec. 4.3: S02 reduction", f"S02 (\\rev{{{f(red['S02'], 1)}\\%}})"),
        ("Fig. 6: S03 L1", f"{f(s03l1[1] - s03l1[0])}~mm higher error for one, S03 L1 ({f(s03l1[1])} against {f(s03l1[0])}~mm"),
        ("Sec. 4.3: S03 L1", f"its error was {f(s03l1[1] - s03l1[0])}~mm higher ({f(s03l1[1])} against {f(s03l1[0])}~mm)"),
        # Fig. 7
        ("Fig. 7: right CL superior capsule", f"errors {f(cl_s.gcz_error_mm.iloc[0])} to {f(cl_s.gczrl_error_mm.iloc[0])}~mm, "
                                              f"{f(cl_s.gcz_error_mm.iloc[1])} and {f(cl_s.gcz_error_mm.iloc[2])}~mm unchanged"),
        ("Fig. 7: CL class", f"CL class error fell from {f(s34.loc['CL', g])} to {f(s34.loc['CL', rl])}~mm"),
        ("Fig. 7: ISL", f"ISL ({f(s34.loc['ISL', g])} to {f(s34.loc['ISL', rl])}~mm)"),
        ("Fig. 7: SSL", f"SSL ({f(s34.loc['SSL', g])} to {f(s34.loc['SSL', rl])}~mm)"),
        ("Fig. 7: ITL tips", f"right {f(itl.loc['right', 'gcz_error_mm'])} to {f(itl.loc['right', 'gczrl_error_mm'])}~mm, left "
                             f"{f(itl.loc['left', 'gcz_error_mm'])} to {f(itl.loc['left', 'gczrl_error_mm'])}~mm; ITL class error "
                             f"{f(s34.loc['ITL', g])} to {f(s34.loc['ITL', rl])}~mm"),
        ("Fig. 7: projection only", f"projection onto the ITL region alone gives {f(stage('projection only (no residual)'))}~mm"),
        ("Fig. 7: inferior CL landmarks", f"moved them {f(cl_big.loc['right', 'correction_mm'], 2)}~mm (right) and "
                                          f"{f(cl_big.loc['left', 'correction_mm'], 2)}~mm (left)"),
        ("Fig. 7: inferior CL errors", f"from {f(cl_big.loc['right', 'gcz_error_mm'], 2)} to {f(cl_big.loc['right', 'gczrl_error_mm'], 2)}~mm "
                                       f"and from {f(cl_big.loc['left', 'gcz_error_mm'], 2)} to {f(cl_big.loc['left', 'gczrl_error_mm'], 2)}~mm"),
        ("Fig. 7: macro errors", f"{f(s34[g].mean())}~mm (GCZ), {f(s34[rl].mean())}~mm (GCZ-RL), and {f(s34['SLD'].mean())}~mm (SLD)"),
        ("Fig. 7: cohort means", f"cohort means of {f(nums['GCZ']['macro_mae'])}, {f(nums['GCZ-RL']['macro_mae'])}, and "
                                 f"{f(nums['SLD']['macro_mae'])}~mm"),
        ("Fig. 7: rank", f"improvement ranks {ords[rank_s03l4]} of the 14 held-out vertebrae"),
        ("Fig. S7: ITL", f"ITL improves from {f(s34.loc['ITL', g], 2)} to {f(s34.loc['ITL', rl], 2)}~mm"),
        # Section 4.4 and Fig. 10
        ("Sec. 4.4: ITL on S03 L4", f"projection alone reduces the ITL error from {f(stage('1 GCZ landmark'))} to "
                                    f"{f(stage('projection only (no residual)'))}~mm and the complete correction to "
                                    f"{f(stage('6 projected = GCZ-RL'))}~mm"),
        ("Fig. 10: stored run", f"GCZ {f(nums['GCZ']['macro_mae'])}~mm, GCZ-RL {f(nums['GCZ-RL']['macro_mae'])}~mm"),
        ("Fig. 10: ten runs", f"GCZ-RL {f(seeds['GCZRL_mm']['min'], 2)}--{f(seeds['GCZRL_mm']['max'], 2)}~mm, mean "
                              f"{f(seeds['GCZRL_mm']['mean'])}~mm; GCZ mean {f(seeds['GCZ_mm']['mean'])}~mm"),
        ("Fig. 10: reduction", f"by {f(seeds['reduction_pct']['min'], 1)}--{f(seeds['reduction_pct']['max'], 1)}\\% (mean "
                               f"{f(seeds['reduction_pct']['mean'], 1)}\\%)"),
        ("Sec. 4.4: seeds", f"GCZ-RL reached {f(seeds['GCZRL_mm']['min'], 2)}--{f(seeds['GCZRL_mm']['max'], 2)}~mm (mean "
                            f"${f(seeds['GCZRL_mm']['mean'], 2)}\\pm{f(seeds['GCZRL_mm']['sd'], 2)}$~mm) and GCZ "
                            f"{f(seeds['GCZ_mm']['min'], 2)}--{f(seeds['GCZ_mm']['max'], 2)}~mm (mean "
                            f"${f(seeds['GCZ_mm']['mean'], 2)}\\pm{f(seeds['GCZ_mm']['sd'], 2)}$~mm)"),
        ("Sec. 4.4: atlas changes", f"Between {int(seeds['atlas_changed']['min'])} and {int(seeds['atlas_changed']['max'])} of the 14"),
        ("Sec. 4.4: projection only", f"{f(abl['GCZ + projection only (no residual)'])}~mm"),
        ("Sec. 4.4: GCZ-RL vs projection only", f"{f(proj['difference_mm'])}"),
        ("Sec. 4.4: interval", f"{f(proj['bootstrap_ci95_mm'][0])}"),
        ("Sec. 4.4: alpha = 1", f"{f(abl['GCZ + bounded residual + projection'])}~mm"),
        ("Table S10: ITL without side constraint", f"{f(rev['no_side_constraint']['MAE_ITL'])}"),
        ("Sec. 4.4: directional agreement", f"{f(rev['interpretability']['pct_cosine_positive_all'], 1)}\\%"),
        # Supplementary S14, S16, S19
        ("S14: ITL footprint in region", f"contains {f(reg.loc['ITL', 'footprint_in_region_pct'], 0)}\\% of its area"),
        ("S14: ITL points inside", f"{f(reg.loc['ITL', 'inside_selection_pct'], 1)}\\% of the ITL points"),
        ("S16: superior above inferior", f"by a median of {f(margin['ALL']['sup_minus_inf_mm_median'], 1)}~mm for ALL and "
                                         f"{f(margin['PLL']['sup_minus_inf_mm_median'], 1)}~mm for PLL (minimum "
                                         f"{f(min(margin['ALL']['sup_minus_inf_mm_min'], margin['PLL']['sup_minus_inf_mm_min']), 1)}~mm)"),
        ("S16: mesh margins", f"within {f(max(abs(v) for d in margin['figure_case_mesh_margin_minus_expert_mm'].values() for v in d.values()), 1)}~mm"),
        ("S16: inferior ALL margin", f"between $z_n={f(lo_i.min(), 2)}$ and ${f(lo_i.max(), 2)}$"),
        ("S16: superior ALL margin", f"the superior between ${f(lo_s.min(), 2)}$ and ${f(lo_s.max(), 2)}$"),
        ("S16: middle half", f"A median of {f(margin['ALL']['region_area_middle_half_pct_median'], 0)}\\% (ALL) and "
                             f"{f(margin['PLL']['region_area_middle_half_pct_median'], 0)}\\% (PLL)"),
        ("S16: separate bands", f"about {f(100 - margin['PLL']['region_area_kept_by_split_bands_pct_median'], 0)}\\% of the PLL region, but "
                                f"{margin['ALL']['gczrl_outside_split_bands'].split('/')[0]} ALL and "
                                f"{margin['PLL']['gczrl_outside_split_bands'].split('/')[0]} PLL landmarks"),
        ("S19: right residual", f"moved the right landmark {f(trace['raw_residual_length_mm']['right'], 2)}~mm towards the tip (error "
                                f"{f(stage('1 GCZ landmark', 'right'), 2)} to {f(stage('2 GCZ + raw residual', 'right'), 2)}~mm)"),
        ("S19: left residual", f"({f(stage('1 GCZ landmark', 'left'), 2)} to {f(stage('2 GCZ + raw residual', 'left'), 2)}~mm)"),
        ("S19: radius", f"radius {f(trace['correction_radius_mm'], 2)}~mm"),
        ("S19: projection only and final", f"reduces the class error to {f(stage('projection only (no residual)'))}~mm, and the complete "
                                           f"correction to {f(stage('6 projected = GCZ-RL'))}~mm"),
        ("S19: residual share", f"the learned residual adds {f(stage('projection only (no residual)') - stage('6 projected = GCZ-RL'))}~mm"),
        ("S19: configuration", f"Extra Trees configuration {trace['regressor']}, $\\alpha_L={trace['alpha']:g}$"),
        # Supplementary S19, independent checks (paper/itl_double_check.py)
        ("S19 check: recomputed errors", "equal those of the evaluation (largest difference below $10^{-13}$~mm)"
         if max((dct[f"{m}_mae_recomputed"] - dct[f"{m}_mae_metrics_file"]).abs().max() for m in (g, rl, "SLD")) < 1e-13 else "(differs)"),
        ("S19 check: tips", f"the expert ITL points lie {f(tip[(g, 'R')]['reference_to_tip_mm'], 2)} (right) and "
                            f"{f(tip[(g, 'L')]['reference_to_tip_mm'], 2)}~mm (left) from the most lateral vertex of their transverse "
                            f"process, and the GCZ-RL landmarks {f(tip[(rl, 'R')]['point_to_tip_mm'], 2)} and {f(tip[(rl, 'L')]['point_to_tip_mm'], 2)}~mm"),
        ("S19 check: GCZ position", f"the GCZ landmarks lie {f(tip[(g, 'R')]['point_to_tip_mm'], 2)} and {f(tip[(g, 'L')]['point_to_tip_mm'], 2)}~mm "
                                    "from it, on surface that faces posteriorly, whereas the surface at both GCZ-RL landmarks faces laterally"
         if all(p["surface_faces"] == ("posterior" if p["method"] == g else "lateral") for p in dc["S03_L4_points"]) else "(differs)"),
        ("S19 check: leakage", f"the {dc['leakage_test_S03']['points_true_run']} S03 landmarks of GCZ and GCZ-RL were identical in both runs and "
                               "equal to the saved predictions"
         if dc["leakage_test_S03"]["max_difference_noisy_vs_true_mm"] == 0 and dc["leakage_test_S03"]["max_difference_recomputed_vs_saved_mm"] < 1e-9
         else "(differs)"),
        ("S19 check: side constraint", f"leaves the ITL errors of all {words[len(dc['S03_ITL_without_side_constraint'])]} S03 vertebrae unchanged"
         if all(v["with_constraint"] == v["without_constraint"] for v in dc["S03_ITL_without_side_constraint"].values()) else "(differs)"),
        ("S19 check: SLD", f"SLD, which shares no step with GCZ-RL, reaches {f(dc['S03_L4']['SLD'])}~mm on this vertebra"),
        ("S19 origin: atlas", f"lie {f(og['S03_L4']['R']['reach_difference_mm'], 2)} (right) and {f(og['S03_L4']['L']['reach_difference_mm'], 2)}~mm "
                              f"(left) closer to the midline than those of S03 L4"),
        ("S19 origin: transferred", f"the transferred points lie {f(og['S03_L4']['R']['transferred_point_medial_of_reference_mm'], 2)} and "
                                    f"{f(og['S03_L4']['L']['transferred_point_medial_of_reference_mm'], 2)}~mm medial of the references"),
        ("S19 origin: recomputed", f"reproduces all {og['points']} saved GCZ ITL points" if og["resnapped_equal_saved"] == og["points"] else "(differs)"),
        ("S19 origin: correlation", f"Over the {og['wide_region_points']} ITL points with a wide selection region, the GCZ error correlates with how much "
                                    f"farther the target's process reaches than the atlas's ($r={f(og['r_error_vs_shortfall_wide'], 2)}$), but not with "
                                    f"the absolute difference ($r={f(og['r_error_vs_abs_difference_wide'], 2)}$)"),
        ("S19 origin: pair", f"the GCZ ITL errors of S01 L1 are {f(og['S01_L1']['R']['gcz_error_mm'], 2)} and {f(og['S01_L1']['L']['gcz_error_mm'], 2)}~mm; "
                             f"with atlas S01 L1, those of S03 L1 are {f(og['S03_L1']['R']['gcz_error_mm'], 2)} and {f(og['S03_L1']['L']['gcz_error_mm'], 2)}~mm"
         if og["S01_L1"]["R"]["atlas"] == "S03_L1" and og["S03_L1"]["R"]["atlas"] == "S01_L1" else "(differs)"),
        ("S19 origin: projection region", f"has {og['S03_L4']['R']['projection_region_vertices_side']} (right) and "
                                          f"{og['S03_L4']['L']['projection_region_vertices_side']} (left) vertices, all within "
                                          f"{f(og['S03_L4']['R']['projection_region_max_to_tip_mm'], 2)} and "
                                          f"{f(og['S03_L4']['L']['projection_region_max_to_tip_mm'], 2)}~mm of the tip"),
        ("S19 check: rank", f"The reduction on S03 L4 ranks {ords[dc['S03_L4_rank_of_ITL_reduction']]} of the {len(dct)} held-out vertebrae; "
                            f"the {words[len(bigger)]} larger ones are {bigger_names}, from {f(bigger[g + '_mae_recomputed'].min(), 2)}--"
                            f"{f(bigger[g + '_mae_recomputed'].max(), 2)}~mm to {f(bigger[rl + '_mae_recomputed'].min(), 2)}--"
                            f"{f(bigger[rl + '_mae_recomputed'].max(), 2)}~mm"),
    ]
    rows = [{"claim": c, "passage": p, "found": p in TEX} for c, p in claims]
    tab = pd.DataFrame(rows)
    tab.to_csv(OUT / "revision" / "text_number_check.csv", index=False)
    bad = tab[~tab.found]
    for r in bad.itertuples():
        print(f"NOT FOUND  {r.claim}:  {r.passage}")
    print(f"{int(tab.found.sum())}/{len(tab)} passages of the revision found verbatim with the regenerated values")
    sys.exit(1 if len(bad) else 0)


if __name__ == "__main__":
    main()
