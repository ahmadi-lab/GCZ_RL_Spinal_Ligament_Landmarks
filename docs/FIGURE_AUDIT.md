# Figure audit

Every figure of the revised manuscript is drawn by `paper/figures_main.py` from saved outputs of
the code release and from the public meshes and annotations (`scripts/fetch_data.py`); no value
is typed in by hand and no anatomy is synthetic. The script also writes
`paper/outputs/figures/figure_audit.csv`, one row per figure with its input files.

Style: white background, Arial, thin axes, light grid; ligament classes in the Okabe-Ito palette;
methods by colour and marker shape (GCZ square, GCZ-RL circle, SLD triangle), validated for colour
vision deficiency with marker shape as the second cue. Quantitative figures are vector PDF and SVG
(text kept as text); figures with 3D anatomy are 600 dpi renders placed in a PDF. All renders use
the same camera conventions (LPS frame; superior, left lateral, and posterior views).

| Figure | Content | Inputs | Decision |
|---|---|---|---|
| Fig. 1 | Pipeline on held-out S03 L4: calibration annotations (S01 L1), calibrated regions, ICP transfer of the atlas, GCZ landmarks, GCZ-RL with expert reference | `P2/outputs/predictions_gcz.csv`, `predictions_gczrl.csv`; raw meshes and annotations of S01, S02, S03 L4 | Replaces the drawn workflow diagram. The ICP transfer is recomputed with the seed of the primary run and reproduces the saved GCZ landmarks of S03 L4 (maximum deviation 6e-14 mm). |
| Fig. 2 | Preprocessing montage | original analysis notebook | Kept unchanged. |
| Fig. 3 | PLL calibration histogram, y_n bands of all classes, candidate regions in three views | raw annotations of S01 L1–L5 and S02 L1–L4; S03 L4 mesh; `gcz.build_thresholds`, `candidate_vertices` | Replaces the drawn calibration figure (which had a smoothed density and a garbled label). |
| Fig. 4 | GCZ landmarks, GCZ-RL landmarks, applied corrections, expert reference on S03 L4; per-class error | `predictions_gcz.csv`, `predictions_gczrl.csv`, `per_vertebra_per_ligament_metrics.csv`; S03 L4 mesh | Replaces the drawn residual-learning diagram. |
| Fig. 5 | Macro MAE by held-out subject and by held-out vertebra | `P2/outputs/table_P2_4_per_subject.csv`, `paper/outputs/primary_per_vertebra_macro_mae.csv` | New. |
| Fig. 6 | GCZ, GCZ-RL, SLD on S03 L4 in three views | `predictions_{gcz,gczrl,sld}.csv`; S03 L4 mesh | Restyled; same points as before (`paper/check_fig6.py`). |
| Fig. 7 | Ablation, coverage and compactness, permutation importance, atlas residual vs error | `paper/outputs/revision/{table_ablation_components,table_region_compactness,feature_importance_permutation,table_atlas_icp}.csv` | Restyled. |
| Fig. 8 | Ten-seed ICP sensitivity (paired GCZ and GCZ-RL per run, stored run, SLD reference) | `paper/outputs/icp_seed_sensitivity/summary.{csv,json}`, `paper/outputs/paper_numbers.json` | Moved from the supplement. Labelled as a computational sensitivity analysis on the same three subjects, not independent validation. |
| Fig. 9 | Development analysis by class (values of Table 9) | `part1_dissertation_57/outputs/table_4_4_heldout_results.csv` | Redrawn as a dot plot without per-bar labels. |
| Fig. S3 | Predicted and true residual vectors on S03 L4 | `paper/outputs/revision/residual_vectors_heldout.csv`; S03 L4 mesh | Restyled. |

The previous Fig. 4 (atlas-to-target prediction workflow) is replaced by panels (c)–(d) of Fig. 1,
and the previous Supplementary Fig. S4 by Fig. 8. The replaced figure files are kept in
`manuscript/figs_previous/`.

Reproduce: run `part2_patient_loso_14/run_part2.py`, `paper/make_paper_numbers.py`,
`paper/revision_analyses.py`, and `paper/icp_seed_sensitivity.py run 0 … 9` + `summary`, then
`python paper/figures_main.py`.
