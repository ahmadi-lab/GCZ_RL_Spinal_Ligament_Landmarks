# Figure audit

The original figures of the manuscript are kept unchanged. New figures were added for the
revision; every new figure is drawn from saved outputs of the code release and from the public
meshes and annotations (`scripts/fetch_data.py`). No value is typed in by hand and no anatomy is
synthetic. `paper/figures_main.py` and `paper/region_verification.py` draw them;
`paper/outputs/figures/figure_audit.csv` lists the input files per figure.

Style of the new figures: white background, Arial, thin axes, light grid; ligament classes in the
Okabe-Ito palette; methods by colour and marker shape (GCZ square, GCZ-RL circle, SLD triangle).
Quantitative figures are vector PDF and SVG; figures with 3D anatomy are 600 dpi renders in a PDF.

## Main text

| Figure | Content | Origin | Inputs |
|---|---|---|---|
| Fig. 1 | GCZ-RL workflow | original (drawn) | — |
| Fig. 2 | Preprocessing montage | original (analysis notebook) | — |
| Fig. 3 | Threshold calibration and candidate regions | original (drawn) | — |
| Fig. 4 | Atlas-to-target prediction workflow | original | — |
| Fig. 5 | Geometry-constrained GCZ-RL correction | original (drawn) | — |
| Fig. 6 | Macro MAE by held-out subject and vertebra | new, `figures_main.py fig5` | `P2/outputs/table_P2_4_per_subject.csv`, `paper/outputs/primary_per_vertebra_macro_mae.csv` |
| Fig. 7 | GCZ, GCZ-RL, SLD on S03 L4 in three views | original (`paper/fig6_qualitative_loso.py`) | saved predictions |
| Fig. 8 | Expert annotations, candidate regions, GCZ-RL landmarks per class on S03 L4 | new, `region_verification.py` | raw meshes and annotations; `gcz.build_thresholds`, `candidate_vertices`; `P2/outputs/predictions_gczrl.csv` |
| Fig. 9 | Component analyses (ablation, compactness, importance, atlas) | new, `figures_main.py fig7` | `paper/outputs/revision/*.csv` |
| Fig. 10 | Ten-seed ICP sensitivity | new, `figures_main.py fig8` | `paper/outputs/icp_seed_sensitivity/summary.{csv,json}`, `paper_numbers.json` |
| Fig. 11 | Development analysis chart | original | values of Table 9 |

## Supplementary material

| Figure | Content | Origin | Inputs |
|---|---|---|---|
| Fig. S1, S2 | Matching diagnostic | original (analysis notebook) | — |
| Fig. S3 | Predicted and true residual vectors on S03 L4 | new, `figures_main.py figS3` | `paper/outputs/revision/residual_vectors_heldout.csv` |
| Fig. S4 | Workflow of Fig. 1 on real data (S03 L4) | new, `figures_main.py fig1` | saved GCZ and GCZ-RL predictions; raw meshes; the ICP transfer reproduces the saved GCZ landmarks (6e-14 mm) |
| Fig. S5 | Correction of Fig. 5 on real data (S03 L4) | new, `figures_main.py fig4` | saved predictions; `per_vertebra_per_ligament_metrics.csv` |

## Candidate-region figure (Fig. 8)

The region figure was checked before it was added (`docs/REGION_VERIFICATION.md`):

- the rendering colours exactly the vertices of each region mask;
- the regions are the masks of the method, unchanged;
- each class column uses one camera for the three rows.

The broad regions of ALL, PLL, and CL are a property of the calibration, which pools the superior
and inferior annotation parts into one band, and are stated as such in the caption.

`figures_main.py` also draws restyled alternatives of Figs. 3, 7, and 11 (`fig3_calibration`,
`fig6_qualitative`, `fig9_development`) that are not used in the manuscript.
