# Part 1 — Dissertation Chapter 4 (57 vertebrae, 20 calibration / 37 held-out)

This folder regenerates the ligament-landmark results reported in Chapter 4 of the
dissertation from the raw vertebral meshes and expert annotations in `../dataset/`.

## 1. How to run

From the package root (`Ligament_Final`), in PowerShell:

```powershell
.\RUN_ALL.ps1 -Part 1            # about 30 min; GCZ geometry read from checkpoints
.\RUN_ALL.ps1 -Part 1 -Fresh     # about 70 min; GCZ geometry recomputed as well
```

or directly:

```powershell
.\.venv\Scripts\python.exe part1_dissertation_57\run_part1.py [--fresh]
```

Everything that is learned (zones, residual targets, correction radii, regressors) is computed
on every run. Only the GCZ geometry (ICP alignment results) can be read from `checkpoints/`;
`--fresh` recomputes it.

## 2. Tables and figures → code → output file

| Dissertation item | What it shows | Produced by (`run_part1.py`) | Output file (`outputs/`) |
|---|---|---|---|
| Table 4.1 | Dataset composition and split | `table_4_1()` | `table_4_1_dataset.csv`, `table_4_1_dataset_per_vertebra.csv` |
| Table 4.2 (= paper Table 3) | Pooled metrics of three matching strategies | `original_notebook/run_original_notebook.py` (original notebook, cells 39 and 41) + `extract_tables.py` | `original_notebook/outputs/tables/table_3_correspondence_strategies__dissertation_table_4_2.csv` |
| Paper Table 1 | Dataset organisation | original notebook, cell 9 | `original_notebook/outputs/tables/table_1_dataset.txt` |
| Paper Table 2 | Bone alignment after ICP (mean NND 1.84 mm) | original notebook, cell 34 | `original_notebook/outputs/tables/table_2_alignment_per_case.csv` |
| Paper Table 4 | Per-ligament MAE, initial vs final matching | original notebook, cell 41 | `original_notebook/outputs/tables/table_4_per_ligament_mae.csv` |
| Paper Table 5 | GCZ vs SLD per ligament | `paper_table_5/run_paper_table_5.py` (original `run_evaluation.py`) | `paper_table_5/outputs/table5_check.csv` |
| — | All paper Tables 1–4 against published values (44/44 identical) | `extract_tables.py` | `original_notebook/outputs/tables/check_against_paper.csv` |
| Table 4.3 | Calibration error of the candidate residual regressors | step 4 (`TABLE_4_3_MODELS`, `gczrl.held_out_records`) | `table_4_3_regressor_comparison.csv`, `table_4_3_full_grid.csv` |
| Table 4.4 | Held-out MAE and Pass@5 per ligament (GCZ, SLD, GCZ-RL) | steps 1, 2, 5, 6 | `table_4_4_heldout_results.csv` |
| Figure 4.10 | Bar chart of Table 4.4 | `figure_4_10()` | `figure_4_10_heldout_results.png/.pdf` |
| Figure 4.11 | GCZ-RL landmarks on a held-out vertebra (S03 L3) | `figure_4_11()` | `figure_4_11_landmarks_S03_L3.png/.csv` |
| Figures 4.1–4.9 | Method illustrations and schematic diagrams | not data figures; the quantities they illustrate are computed by `code/ligament/gcz.py` (zones, ICP transfer, GCZ point) and `code/ligament/gczrl.py` (features, residual, clipping, projection) | — |
| — | Regenerated vs dissertation values; regenerated vs archived predictions | step 6 | `verification_part1.json` |

Per-vertebra results: `per_vertebra_per_ligament_metrics.csv`. Predicted coordinates:
`predictions_gcz_heldout37.csv`, `predictions_gczrl_extra_trees_heldout37.csv`.
Trained model: `models/gczrl_extra_trees_final.joblib`.

## 3. Data and split

* 57 vertebral surface meshes from five subjects (S01–S05); `../dataset/split_manifest_57.csv`
  fixes the split.
* Calibration (20): S01 L1–L5 and S02 C3–T10.
* Held-out (37): S02 T11–L4 (6), S03 C3–L5 (21), S04 C3–C7 (5), S05 C3–C7 (5).
* **All 57 vertebrae carry expert annotations.** S01–S03 (47 vertebrae) are in
  `../dataset/raw/VerSe_Ligaments_GT`; S04–S05 (10 cervical vertebrae, Bc11701/Bc11702) are in
  `../dataset/raw/Cervical_Ligaments_GT` (public SpineToolkit `Cervical_Ligaments_GT`, 15 markup
  files per vertebra, meshes identical to the package meshes).
* The original evaluation read annotations from the VerSe folder only (`BASE_CERV = BASE` in the
  original notebook), so Table 4.4 scores **27 of the 37 held-out vertebrae** (S02: 6, S03: 21);
  the 10 held-out cervical vertebrae are predicted but not scored. `run_part1.py` keeps this so
  that Table 4.4 matches the dissertation.
* The split is at the level of the vertebra. Subject S02 contributes vertebrae to both
  calibration (15) and held-out (6) sets. Part 2 provides the patient-level evaluation.

## 4. Method (single final version: GCZ-RL with Extra Trees)

1. **GCZ (geometry only).** Zones per ligament from the calibration annotations
   (AABB-normalised coordinates, mean ± 2 SD); atlas = the five S01 lumbar calibration vertebrae;
   best atlas by quick ICP, full ICP, transfer of the atlas landmarks, snapping into the zone,
   10 points per ligament (`code/ligament/gcz.py`).
2. **Training data (calibration vertebrae only).** GCZ is run on each calibration vertebra
   (an atlas vertebra is never used to predict itself); every GCZ point is paired with an expert
   landmark (Hungarian assignment); the difference is the residual target. Bilateral ligaments
   are mirrored across the midline. Result: **1683 samples × 36 geometric features** from the
   20 calibration vertebrae.
3. **Residual regressor.** One Extra Trees regressor per ligament class
   (`FINAL_CONFIGURATION` in `code/ligament/gczrl.py`):

   | Ligament | Extra Trees configuration | Blend weight α |
   |---|---|---|
   | ALL | ET-2: 520 trees, unlimited depth, all features | 1.20 |
   | PLL | ET-2 | 1.20 |
   | LF | ET-1: 220 trees, depth 18, 85 % features | 1.00 |
   | ISL | ET-4: 520 trees, depth 18, 85 % features, bootstrap | 1.20 |
   | SSL | ET-1 | 1.00 |
   | ITL | ET-4 | 1.00 |
   | CL | ET-2 | 1.10 |

4. **Correction.** Predicted residual clipped to the correction radius τ (mean + 2 SD of the
   calibration residual magnitudes of that class), scaled by α, and projected onto the nearest
   vertex inside the ligament zone.
5. The models are trained once on the 20 calibration vertebrae and applied once to the 37
   held-out vertebrae.

## 5. How the final configuration was fixed (read before citing the 20/37 numbers)

* Table 4.3 compares Ridge, KNN, Random Forest, Gradient Boosting and Extra Trees by
  leave-one-vertebra-out error on the 20 calibration vertebrae, for blend weights 0.70–1.20.
  Each column is the best configuration of that family (Extra Trees: best of ET-1…ET-4;
  Gradient Boosting: best of the two boosting implementations).
* For ALL, PLL, ISL, ITL and CL the final entry is the calibration-best Extra Trees configuration
  and α of Table 4.3.
* For LF and SSL the final entry (ET-1, α = 1.00) is not the calibration-best entry of Table 4.3.
  It was retained after the held-out results of the candidate configurations had been seen
  during development. These two choices therefore used held-out information, and the 20/37
  numbers should be read as optimistic. Part 2 makes every choice on training subjects only and
  is the unbiased estimate.
* Code check: in Table 4.3 the Random-Forest value for LF is 3.329 mm in the code output and
  3.359 mm in the dissertation; with the code value Random Forest, not Extra Trees, has the
  lowest calibration error for LF. The original selection script of the study
  (`model_comparison_per_ligament.csv` of the archived results) gives the same 3.329 mm and selects
  Random Forest for LF. Gradient Boosting for LF is 3.545 mm in the code and 3.544 mm in the
  dissertation (rounding). All other cells agree with the dissertation.
* Computer dependence of Table 4.3: the Gradient Boosting column is not bit-identical across
  computers. On the computer used for the current outputs it differs from the values above by at
  most 0.029 mm (ALL 1.627, PLL 1.583, LF 3.543, ISL 1.919, CL 5.157 mm; 35 of 315 grid rows,
  all Gradient Boosting). No family ranking changes, and the Extra Trees, Ridge, KNN and Random
  Forest columns, Table 4.4 and all predictions are reproduced exactly.

## 6. Notes on individual tables

* **Table 4.4** is reproduced exactly from the code (all MAE and Pass@5 values, all ligaments).
  The metrics use one-to-one Hungarian matching between predicted and expert points
  (`code/ligament/metrics.py`). MAE = mean over vertebrae and then over the seven ligaments;
  Pass@5 = matched pairs within 5 mm, pooled within a ligament and averaged over ligaments.
  (The dissertation text describes this table as using the anatomy-constrained matching of
  Section 4.3.3.1; the reproduced values are those of plain Hungarian matching.)
* **Table 4.4, GCZ column** comes from a GCZ run on the STL surfaces (step 1). GCZ-RL corrects a
  second GCZ run on the VTK surfaces (step 2; same algorithm, mean error 3.704 mm vs 3.708 mm).
  Both runs are regenerated exactly.
* **Table 4.4, SLD column** is computed from the SLD predictions saved for the dissertation
  (`../dataset/dissertation_archived_predictions/sld_heldout37.csv`). That SLD run used the
  surface coordinates without conversion to the coordinate frame SLD expects; the corrected SLD
  execution is reported in Part 2.
* **Table 4.2 and paper Tables 1–4** come from the original GCZ notebook
  `original_notebook/GCZ_journal_FINAL_57_20_37_eval9_ORIGINAL.ipynb` (unchanged; SHA-256
  `691f3ee2…`). `run_original_notebook.py` executes it on the packaged data, changing only the
  data/output paths and not opening a browser. All 37 code cells run without error; 32 give text
  output identical to the saved notebook and the other 5 differ only in file paths, library
  versions, or a figure saved as HTML instead of PNG (no `kaleido`). `extract_tables.py` reproduces
  all 44 published values of paper Tables 1–4. These tables use the 9 primary lumbar held-out
  vertebrae (Verse561 L1–L4, Verse581 L1–L5) and the notebook's own GCZ run (1126 Hungarian
  pairs), not the GCZ of Table 4.4 (27 vertebrae). The labels inside the notebook
  (v4, v5, … v8d) name successive *matching rules* tried in that notebook, not versions of
  GCZ-RL. Table 4.2 rows: "Pure Euclidean Hungarian" = cell 36 strategy A; "Initial anatomy-aware
  matching" = cell 34 (rule v5); "Final anatomy-constrained GCZ" = cell 41 (rule v8d).
* **Table 4.1**: all 57 vertebrae carry annotations, as the dissertation states. The output also
  lists how many were used in the evaluation (20 calibration + 27 held-out; Section 3).

* **Paper Table 5** is regenerated by `paper_table_5/run_paper_table_5.py`, which runs the original
  evaluation script with packaged paths; every per-ligament value matches the paper, both with the
  archived GCZ predictions and with GCZ regenerated by `run_part1.py`. It averages per
  vertebra–ligament case over the cases that have an annotation of that ligament (SSL n = 14,
  ITL n = 21), so its values differ from Table 4.4. The overall SLD Pass@5 printed in the paper (83.3 %) differs
  from the mean of its seven per-ligament values (84.3 %).
* `table_4_2_source/` holds a second copy of the same notebook (identical code and outputs,
  different metadata), kept for reference.

## 7. Data-leakage controls

* Held-out annotations are read only in the evaluation step: during training the label function
  returns nothing for held-out vertebrae (`gt_calib` in `run_part1.py`).
* Zones, atlas, benchmark prior, residual targets, correction radii and regressors are computed
  from the 20 calibration vertebrae only.
* Regenerated GCZ and GCZ-RL predictions are compared point by point with the predictions saved
  for the dissertation (`verification_part1.json`, `regenerated_vs_archived_predictions`).
* Limitation of this protocol: the split is by vertebra, so S02 appears on both sides
  (Section 3), and two configuration entries used held-out information (Section 5). Both are
  removed in Part 2.

## 8. Limitations

* The experimental data are the vertebrae available for this study: three annotated subjects
  (S01–S03) from one public CT dataset (VerSe) with one set of expert annotations. No second
  annotator, no pathology labels, and no subject-level confidence intervals.
* Landmark accuracy is not ligament-reconstruction accuracy; the effect on finite element
  predictions has not been measured.
