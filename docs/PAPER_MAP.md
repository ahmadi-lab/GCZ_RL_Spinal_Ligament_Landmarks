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
| Fig. 1 | Pipeline on held-out S03 L4 (calibration, regions, ICP transfer, GCZ, GCZ-RL) | `paper/figures_main.py` `fig1` | `paper/outputs/figures/fig1_pipeline.{pdf,png}` |
| Fig. 2 | Preprocessing; panel D midplane positions | `NB` cells 21–24 | figure only |
| Fig. 3 | Calibration bands and candidate regions (fold holding out S03) | `paper/figures_main.py` `fig3` | `paper/outputs/figures/fig3_calibration.{pdf,png}` |
| Fig. 4 | GCZ initialization, GCZ-RL correction, expert reference on S03 L4 | `paper/figures_main.py` `fig4` | `paper/outputs/figures/fig4_residual_correction.{pdf,png}` |
| Fig. 5 | Primary evaluation by held-out subject and vertebra | `paper/figures_main.py` `fig5` | `paper/outputs/figures/fig5_patient_level.{pdf,svg,png}` |
| Fig. 6 | Qualitative, three orthogonal views | `paper/figures_main.py` `fig6` | `paper/outputs/figures/fig6_qualitative.{pdf,png}` |
| Fig. 6 caption | Errors on S03 L4 (3.263 / 1.972 / 2.771 mm), rank 5 of 14; check that the drawn points are the saved predictions | `paper/check_fig6.py` | `paper/outputs/fig6_check_S03_L4.csv`, `fig6_vertebra_ranking.csv` |
| Fig. 7 | Component analyses (ablation, compactness, importance, atlas) | `paper/figures_main.py` `fig7` | `paper/outputs/figures/fig7_components.{pdf,svg,png}` |
| Fig. 8 | Ten-seed ICP sensitivity | `paper/figures_main.py` `fig8` | `paper/outputs/figures/fig8_icp_seeds.{pdf,svg,png}` |
| Fig. 9 | Development analysis chart (values of Table 9) | `paper/figures_main.py` `fig9` | `paper/outputs/figures/fig9_development.{pdf,svg,png}` |
| Graphical abstract | — | `paper/graphical_abstract.py` | `paper/outputs/graphical_abstract.{png,pdf}` |

Table rows of Tables 6–8 and S10–S15 are written by `paper/revision_tables.py` into
`manuscript/tables/` and read by the LaTeX source.

### Numbers in the text

| Section | Numbers | Source |
|---|---|---|
| Sec. 4.3 | statistics (Wilcoxon, cluster bootstrap), 13/14 vertebrae, ITL errors | `paper/outputs/paper_numbers.json` |
| Sec. 4.3 | leakage test (424 points; 412 with recomputed ICP) | `P2/outputs/table_P2_5_leakage_test.json`, `REV/leakage_test_recomputed_icp_this_environment.json` |
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
| Table S16 (ICP seed sensitivity; Sec. 4.4, Fig. 8, Limitations) | `paper/icp_seed_sensitivity.py run 0 … 9`, then `summary` | `paper/outputs/icp_seed_sensitivity/{offset_*.json,summary.csv,summary.json,tab_icp_seeds.tex}` |

Every figure is listed with its input files in `docs/FIGURE_AUDIT.md` and
`paper/outputs/figures/figure_audit.csv` (written by `paper/figures_main.py`).
