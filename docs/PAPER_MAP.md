# Paper item → script → output

Numbering follows the revised journal manuscript. Abbreviations: `P1` =
`part1_dissertation_57`, `P2` = `part2_patient_loso_14`, `NB` =
`P1/original_notebook` (the original analysis notebook, executed unchanged by
`run_original_notebook.py`), `REV` = `paper/outputs/revision`.

## Main text

| Item | Content | Script | Output |
|---|---|---|---|
| Table 1 | Dataset and study design | `P1/run_part1.py`; cohort in `P2/run_part2.py` | `P1/outputs/table_4_1_dataset.csv`, `dataset/split_manifest_57.csv` |
| Table 2 | Alignment quality | `paper/make_paper_numbers.py` (primary); `NB` cell 34 (development) | `paper/outputs/paper_numbers.json` → `atlas_residual_mm`; `NB/outputs/tables/table_2_alignment_per_case.csv` |
| Table 3 | Calibration regressor comparison | `paper/table3_regressor_comparison.py` | `paper/outputs/table3_rows_fixed_ctx.csv` |
| Table 4 | Primary LOSO, per ligament | `P2/run_part2.py`, `paper/make_paper_numbers.py` | `P2/outputs/table_P2_3_per_ligament.csv` (MAE); `paper/outputs/table4_primary_per_ligament.csv` (Pass@5 averaged over cases) |
| Table 5 | Primary LOSO, per subject | `P2/run_part2.py` | `P2/outputs/table_P2_4_per_subject.csv`, `table_P2_1_folds.csv` |
| Table 6 | Primary comparison, matched-pair and matching-independent metrics | `paper/make_paper_numbers.py`, `paper/revision_tables.py` | `paper/outputs/paper_numbers.json` |
| Table 7 | Component ablation | `paper/revision_analyses.py ablation` | `REV/table_ablation_components.csv`, `REV/ablation_per_case.csv` |
| Table 8 | Candidate-region coverage and compactness | `paper/revision_analyses.py compactness` | `REV/table_region_compactness.csv`, `REV/region_compactness_per_vertebra.csv` |
| Table 9 | Development analysis | `P1/run_part1.py` | `P1/outputs/table_4_4_heldout_results.csv` |
| Fig. 1 | GCZ-RL workflow (original drawing) | — (drawn) | manuscript source `figs/fig_gcz_pipeline_benchmark_style4.png` |
| Fig. 2 | Preprocessing; panel D midplane positions | `NB` cells 21–24 | figure only |
| Fig. 3 | Threshold calibration and candidate regions (original drawing) | — (drawn) | `figs/fig_threshold_candidate_combined1.png` |
| Fig. 4 | Atlas-to-target prediction workflow (original) | — | `figs/fig_prediction_pipeline.png` |
| Fig. 5 | Geometry-constrained correction (original drawing) | — (drawn) | `figs/gcz_residual_learning_framework.png` |
| Fig. 6 | Primary evaluation by held-out subject and vertebra | `paper/figures_main.py` `fig5` | `paper/outputs/figures/fig5_patient_level.{pdf,svg,png}` |
| Fig. 7 | GCZ and GCZ-RL on S03 L4 with magnified close-ups of the applied correction; per-class errors of GCZ, GCZ-RL, SLD | `paper/figures_anatomy.py fig_qualitative_insets` | `paper/outputs/figures/fig_qualitative_insets.{pdf,png}`, checks in `anatomy_checks.json` |
| Fig. 7 and Fig. S7 captions | Errors on S03 L4 (3.263 / 1.972 / 2.771 mm), rank 5 of 14; check that the drawn points are the saved predictions | `paper/check_fig6.py` | `paper/outputs/fig6_check_S03_L4.csv`, `fig6_vertebra_ranking.csv` |
| Fig. 8 | Expert annotations, calibrated search regions, GCZ-RL landmarks; ALL and PLL margin close-ups (S03 L4) | `paper/figures_anatomy.py fig_regions_margins` | `paper/outputs/figures/fig_regions_margins.{pdf,png}` (the earlier 3 × 7 grid stays in `paper/outputs/region_verification/`) |
| Fig. 9 | Component analyses (ablation, compactness, importance, atlas) | `paper/figures_main.py` `fig7` | `paper/outputs/figures/fig7_components.{pdf,svg,png}` |
| Fig. 10 | Ten-seed ICP sensitivity | `paper/figures_main.py` `fig8` | `paper/outputs/figures/fig8_icp_seeds.{pdf,svg,png}` |
| Fig. 11 | Development analysis chart (original) | values of Table 9 | `figs/fig_literature_comparison_clean2.png` |
| Graphical abstract | — | `paper/graphical_abstract.py` | `paper/outputs/graphical_abstract.{png,pdf}` |

Table rows of Tables 6–8 and S10–S14 are written by `paper/revision_tables.py` into
`manuscript/tables/`; those of Tables S15–S18 by `reproducibility_check.py`, `icp_seed_sensitivity.py summary`,
`region_verification.py`, and `itl_trace.py`, copied by `paper/build_manuscript_assets.py`. The LaTeX source
reads them. `paper/check_text_numbers.py` checks the numbers that the revision added to the captions and
to Sections 4.3–4.4 and S13–S19 against the regenerated outputs.

### Numbers in the text

| Section | Numbers | Source |
|---|---|---|
| Sec. 4.3 | statistics (Wilcoxon, cluster bootstrap), 13/14 vertebrae, ITL errors | `paper/outputs/paper_numbers.json` |
| Sec. 4.3 | leakage test (424 points; 412 with recomputed ICP) | `P2/outputs/table_P2_5_leakage_test.json` (archived original run, 424), `P2/outputs/table_P2_5_leakage_test_recomputed.json` (written by `run_leakage_test.py` in the current environment; 412 here) |
| Sec. 4.3 | selected α and models | `P2/outputs/table_P2_6_selected_models.csv` |
| Sec. 4.4 | ablation, GCZ-RL vs GCZ with projection, side-constraint re-run, compactness, interpretability, atlas | `REV/revision_numbers.json` |
| Sec. 4.5 | GCZ input of GCZ-RL (3.704 mm, 80.8 %) | `paper_numbers.json` → `dev_GCZ_input_to_GCZRL` |

`paper/make_paper_numbers.py` compares every value listed in `paper/expected_values.json`
with the regenerated value, and exits with an error if any differs.
`paper/reproducibility_check.py` compares 68 main claims with the regenerated outputs
(Supplementary Table S15).

## Supplementary material

| Item | Script | Output |
|---|---|---|
| Table S1 (features) | `code/ligament/gczrl.py` `point_features`, `_zone_features` | — |
| Tables S2–S3 | `paper/make_paper_numbers.py` | `paper/outputs/table_S2_S3_point_accounting_and_matching_independent.csv` |
| Table S4 | `paper/make_paper_numbers.py` | `paper/outputs/table_S4_candidate_region_feasibility.csv` |
| Table S5 (settings) | constants in `code/ligament/gcz.py`, `gczrl.py`, `sld_official.py` | — |
| Table S6 (regressor correction) | `paper/table3_regressor_comparison.py` | `paper/outputs/table3_{as_published,rows_fixed,rows_fixed_ctx}.csv` |
| Tables S7–S8, Figs. S1–S2 (matching diagnostic) | `NB` cells 36–41 | `NB/outputs/tables/table_3_*.csv`, `table_4_per_ligament_mae.csv` |
| Table S9 (selected models) | `P2/run_part2.py` | `P2/outputs/table_P2_6_selected_models.csv` |
| Table S10 (ablation by class, side-constraint re-run) | `paper/revision_analyses.py ablation no_side_constraint` | `REV/table_ablation_components.csv`, `REV/table_ablation_no_side_constraint.csv` |
| Table S11 (compactness by subject) | `paper/revision_analyses.py compactness` | `REV/region_compactness_per_vertebra.csv` |
| Tables S12–S13, Fig. S3 (residual interpretability) | `paper/revision_analyses.py interpretability`, `paper/figures_main.py` `figS3` | `REV/table_residual_agreement.csv`, `REV/feature_importance_*.csv`, `REV/residual_vectors_heldout.csv`, `paper/outputs/figures/figS3_residual_vectors.{pdf,png}` |
| Table S14 (atlas and ICP) | `paper/revision_analyses.py atlas` | `REV/table_atlas_icp.csv` |
| Table S15 (reproducibility) | `paper/reproducibility_check.py` | `REV/reproducibility_table.{csv,tex}` |
| Table S16 (ICP seed sensitivity; Sec. 4.4, Fig. 10, Limitations) | `paper/icp_seed_sensitivity.py run 0 … 9`, then `summary` | `paper/outputs/icp_seed_sensitivity/{offset_*.json,summary.csv,summary.json,tab_icp_seeds.tex}` |

| Table S17 (region verification) | `paper/region_verification.py` | `paper/outputs/region_verification/{region_verification_by_class,region_verification_per_vertebra,bands_fold_S03}.csv`, `tab_region_verification.tex` |
| Figs. S4–S5 (method on real data) | `paper/figures_main.py` `fig1`, `fig4` | `paper/outputs/figures/fig1_pipeline.*`, `fig4_residual_correction.*` |
| Section S16, Fig. S6 (ALL and PLL margins) | `paper/margin_verification.py` | `paper/outputs/margin_verification/{margin_check,margin_zn_by_vertebra}.csv`, `margin_summary.json`, `fig_margins_ALL_PLL.{pdf,png}` |
| Section S17, Fig. S7 (previous Fig. 7, unchanged) | `paper/fig6_qualitative_loso.py` | `paper/outputs/fig6_qualitative_loso_S03_L4.png` |
| Section S18, Fig. S8 (calibrated GCZ search regions) | `paper/figures_anatomy.py` | `paper/outputs/figures/fig_search_regions.{pdf,png}` |
| Section S19, Table S18 (step-by-step ITL correction on S03 L4) | `paper/itl_trace.py` | `paper/outputs/revision/itl_trace_S03_L4.{csv,json}`, `tab_itl_trace.tex` |

Every figure is listed with its input files in `docs/FIGURE_AUDIT.md` and
`paper/outputs/figures/figure_audit.csv` (written by `paper/figures_main.py`); the candidate
regions are checked against the expert annotations in `docs/REGION_VERIFICATION.md`.
