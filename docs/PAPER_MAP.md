# Paper item → script → output

Numbering follows the Medical Image Analysis manuscript. Abbreviations: `P1` =
`part1_dissertation_57`, `P2` = `part2_patient_loso_14`, `NB` =
`P1/original_notebook` (the original analysis notebook, executed unchanged by
`run_original_notebook.py`).

## Main text

| Item | Content | Script | Output |
|---|---|---|---|
| Table 1 | Dataset and study design | `P1/run_part1.py`; cohort in `P2/run_part2.py` | `P1/outputs/table_4_1_dataset.csv`, `dataset/split_manifest_57.csv` |
| Table 2 | Alignment quality | `paper/make_paper_numbers.py` (primary); `NB` cell 34 (development) | `paper/outputs/paper_numbers.json` → `atlas_residual_mm`; `NB/outputs/tables/table_2_alignment_per_case.csv` |
| Table 3 | Calibration regressor comparison | `paper/table3_regressor_comparison.py` | `paper/outputs/table3_rows_fixed_ctx.csv` |
| Table 4 | Primary LOSO, per ligament | `P2/run_part2.py`, `paper/make_paper_numbers.py` | `P2/outputs/table_P2_3_per_ligament.csv` (MAE); `paper/outputs/table4_primary_per_ligament.csv` (Pass@5 averaged over cases) |
| Table 5 | Primary LOSO, per subject | `P2/run_part2.py` | `P2/outputs/table_P2_4_per_subject.csv`, `table_P2_1_folds.csv` |
| Table 6 | Development analysis | `P1/run_part1.py` | `P1/outputs/table_4_4_heldout_results.csv` |
| Fig. 1, 3–5 | Method schematics | — (drawn) | manuscript source `figs/` |
| Fig. 2 | Preprocessing; panel D midplane positions | `NB` cells 21–24 | figure only |
| Fig. 6 | Qualitative, three orthogonal views | `paper/fig6_qualitative_loso.py` | `paper/outputs/fig6_qualitative_loso_S03_L4.png` |
| Fig. 7 | Development analysis chart | values of Table 6 | `P1/outputs/figure_4_10_heldout_results.png` (same data) |
| Graphical abstract | — | `paper/graphical_abstract.py` | `paper/outputs/graphical_abstract.{png,pdf}` |

### Numbers in the text

| Section | Numbers | Source |
|---|---|---|
| Sec. 4.3 | statistics (Wilcoxon, cluster bootstrap), 13/14 vertebrae, ITL errors | `paper/outputs/paper_numbers.json` |
| Sec. 4.3 | leakage test | `P2/outputs/table_P2_5_leakage_test.json` |
| Sec. 4.3 | selected α and models | `P2/outputs/table_P2_6_selected_models.csv` |
| Sec. 4.4 | GCZ input of GCZ-RL (3.704 mm, 80.8 %) | `paper_numbers.json` → `dev_GCZ_input_to_GCZRL` |

`paper/make_paper_numbers.py` compares every value listed in `paper/expected_values.json`
with the regenerated value, and exits with an error if any differs.

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
