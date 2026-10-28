# Figure audit

Figure-by-figure audit of the revised manuscript (2026-10-28). The original figures of the
manuscript are kept unchanged as images; new figures were added for the revision. Every new figure
is drawn from saved outputs of the code release and from the public meshes and annotations
(`scripts/fetch_data.py`). No value is typed in by hand, no anatomy is synthetic, and no mask or
prediction was changed for a figure. `paper/outputs/figures/figure_audit.csv` lists the input files
per figure, and `paper/outputs/figures/anatomy_checks.json` holds the automatic checks of the
anatomical figures. Before-and-after previews of every changed figure, and the compiled pages
before and after, are in `paper/outputs/figure_audit/` and in the PDF report
`GCZ_RL_Figure_Audit.pdf` of the submission package.

## Numbering: compiled figure, source file, script

The compiled numbers differ from several file names, which are historical. Edit by file, not by
number.

| Compiled | File in the manuscript source | Script | Origin |
|---|---|---|---|
| Fig. 1 | `figs/fig_gcz_pipeline_benchmark_style4.png` | — (drawn) | original |
| Fig. 2 | `figs/fig_preprocessing_montage.png` | analysis notebook | original |
| Fig. 3 | `figs/fig_threshold_candidate_combined1.png` | — (drawn) | original |
| Fig. 4 | `figs/fig_prediction_pipeline.png` | — | original |
| Fig. 5 | `figs/gcz_residual_learning_framework.png` | — (drawn) | original |
| Fig. 6 | `figs/fig5_patient_level.pdf` | `paper/figures_main.py fig5` | new |
| Fig. 7 | `figs/fig_qualitative_insets.pdf` | `paper/figures_anatomy.py fig_qualitative_insets` | new (2026-10-28) |
| Fig. 8 | `figs/fig_regions_margins.pdf` | `paper/figures_anatomy.py fig_regions_margins` | new (2026-10-28) |
| Fig. 9 | `figs/fig7_components.pdf` | `paper/figures_main.py fig7` (alternatives: `fig7_alternatives`) | new |
| Fig. 10 | `figs/fig8_icp_seeds.pdf` | `paper/figures_main.py fig8` | new |
| Fig. 11 | `figs/fig11_development_results.pdf` | `paper/fig11_development.py` (Table 9 values) | redrawn from code in the design of the original; the original image `fig_literature_comparison_clean2.png` is kept in `figs/` |
| Fig. S1, S2 | `fig_matching_versions_clean.png` (output of the original notebook, `run_original_notebook.py`), `figure4_cropped_labeled1.png` (drawn) | analysis notebook / drawn | original |
| Fig. S3 | `figS3_residual_vectors.png` | `paper/figures_main.py figS3` | new |

Every generated file is copied into `manuscript/` by `paper/build_manuscript_assets.py`, which fails if a figure or table of `main.tex` or `supplementary.tex` has no source; only the drawn schematics of Figs. 1–5 and S2 are kept as files. `run_all.sh` runs every step; `PDF=1` also compiles the clean and highlighted manuscript and the supplement.
| Fig. S4 | `fig1_pipeline.pdf` | `paper/figures_main.py fig1` | new |
| Fig. S5 | `fig4_residual_correction.pdf` | `paper/figures_main.py fig4` | new |
| Fig. S6 | `fig_margins_ALL_PLL.pdf` | `paper/margin_verification.py` | new |
| Fig. S7 | `fig_qualitative_loso_S03_L4.png` | `paper/fig6_qualitative_loso.py` | original (was Fig. 7) |

`figures_main.py` also draws restyled alternatives (`fig3_calibration`, `fig6_qualitative`,
`fig9_development`) that are not used in the manuscript; `fig6_qualitative.png` in the figure
package is one of these and is **not** the compiled Fig. 7.

## Encoding used in all new figures

- Methods: GCZ blue square, GCZ-RL red circle, SLD green triangle; paired GCZ and GCZ-RL values
  joined by a gray segment.
- Anatomy (Figs. 7, 8): expert attachment annotation black sphere; GCZ landmark blue cube; GCZ-RL
  landmark red sphere; applied correction dark segment; calibrated search region (GCZ mask) gray
  surface tint. The three representations (annotations, search regions, final landmarks) never
  share a colour or a mark.
- Markers in the anatomical figures are glyphs of fixed size in millimetres; every view uses a
  parallel projection with an explicit camera, and a close-up is a crop of the same render as its
  overview, so the frames on the overview are exact.
- Arial, 7–8 pt at print size; line widths 0.8–1.2 pt; legends outside the data.

## Figure by figure

| Figure | Data source | Orientation and coordinates | Landmarks against saved predictions | Print-size check | Change in this audit |
|---|---|---|---|---|---|
| Fig. 1 | drawn | — | — | readable at text width | none (original); real-data version is Fig. S4 |
| Fig. 2 | notebook | — | — | readable | none (original) |
| Fig. 3 | drawn | — | — | readable | caption: the drawn zones are search regions, not attachment footprints; points to Fig. 8 |
| Fig. 4 | original | — | — | readable | none (original) |
| Fig. 5 | drawn | — | — | readable | none (original); real-data version is Fig. S5 |
| Fig. 6 | `table_P2_4_per_subject.csv`, `primary_per_vertebra_macro_mae.csv`, `per_vertebra_per_ligament_metrics.csv` | — | exact values from the per-case errors; assert they round to Table 5 and the per-vertebra file | readable at text width | (a) paired GCZ/GCZ-RL per subject with the reduction (38.8, 10.0, 15.7 %); (b) GCZ minus GCZ-RL per vertebra, exception S03 L1 marked; SLD only in the tables |
| Fig. 7 | `predictions_gcz.csv`, `predictions_gczrl.csv`, `per_vertebra_per_ligament_metrics.csv`, S03 L4 mesh | one posterior camera (patient's right on the right, LPS); close-ups b–e are crops of the overview render, so each frame is exactly its close-up | drawn = saved rows (0 mm); class errors recomputed = metrics file; `fig7_landmark_verification_S03_L4.csv` lists every landmark, matched reference, error, region membership, visibility | CL_S_R references 10/10 visible; ITL tips 2/2 | (b) right facet joint (superior capsule); (c) ISL and SSL; (d, e) full transverse processes (ITL 11.508 to 0.590 and 4.682 to 0.765 mm); posteroinferior view removed; values in the caption |
| Fig. 8 | `predictions_gczrl.csv`, annotation files, masks from `gcz.build_thresholds` and `candidate_vertices` (fold holding out S03) | ALL anterior, PLL posterior view of the body; small regions end-on | the dark-gray bands were the actual masks: they appear only in Part B, unchanged (gray vertices = mask vertices, asserted); Part A footprints = surface within 2.5 mm of the expert points | close-ups at one scale (46 × 21 mm), two frames per vertebra | Part A: expert footprints (superior, inferior, other), expert points, GCZ-RL; Part B: actual masks |
| Fig. 9 | `paper/outputs/revision/*.csv` | — | values plotted from the tables | legends outside the data; axis breaks marked; headings aligned | (a) header and spacing, caption on α_L = 1; (b) labels right of points; (c) fold ticks, mean, range; (d) white edges for coincident markers |
| Fig. 10 | `icp_seed_sensitivity/summary.{csv,json}`, `paper_numbers.json` | — | paired values per run unchanged; GCZ-RL < GCZ asserted for all ten runs | full text width (the numbers beside the old column figure were the line numbers of the other column) | stored run apart (gap and divider), filled markers for the ten runs, SLD value in the legend |
| Fig. 11 | `part1_dissertation_57/outputs/table_4_4_heldout_results.csv` | — | values plotted from the table | vector, readable at 0.85 of the text width | drawn from code; the bar heights of the original image did not match its printed values exactly (e.g. GCZ Pass@5 for ALL: bar about 89 %, label 86.7 %) |
| Fig. S3–S5 | see above | as in the main figures | S4: ICP transfer reproduces the saved GCZ landmarks (6e-14 mm) | checked | S5: wider gap before panel (c) |
| Fig. S6 | `_S`/`_I` annotation files, masks, saved predictions | superior above inferior in 28/28 cases | — | checked | none |
| Fig. S7 | saved predictions | — | `paper/check_fig6.py` | readable | moved from the main text unchanged; caption value 8.09 mm (the previous 8.10 came from rounding 8.095) |

## Values shown in the figures

All values printed in the figures are computed in the plotting scripts from the saved outputs and
equal the manuscript tables: S03 L4 macro MAE 3.263 / 1.972 / 2.771 mm (GCZ / GCZ-RL / SLD); class
errors on S03 L4 from `per_vertebra_per_ligament_metrics.csv`; ablation 2.879, 2.310, 2.698, 2.698,
2.132, 2.150 mm; the ten-seed values of Table S16. One rounding was corrected: the ITL error of
GCZ on S03 L4 is 8.0948 mm and is now given as 8.09 mm (previously 8.10).

## Candidate regions

The broad regions of ALL, PLL, and CL are a property of the calibration and are stated as such in
the captions; nothing was shrunk or moved (`docs/REGION_VERIFICATION.md`, Supplementary S14, S16).
