# Part 2 — Patient-level validation (leave-one-subject-out, 14 lumbar vertebrae)

Part 1 reproduces the dissertation, whose split is by vertebra. This part answers the stricter
question: how does the same method perform on a subject it has never seen?

## 1. How to run

From the package root (`Ligament_Final`):

```powershell
.\RUN_ALL.ps1 -Part 2            # about 7 min; GCZ geometry and SLD output read from checkpoints
.\RUN_ALL.ps1 -Part 2 -Fresh     # about 10 min; GCZ geometry and SLD recomputed as well
```

or directly:

```powershell
.\.venv\Scripts\python.exe part2_patient_loso_14\run_part2.py [--fresh]
.\.venv\Scripts\python.exe part2_patient_loso_14\run_leakage_test.py
```

Model training and model selection are executed on every run. Checkpoints contain only the GCZ
geometry (ICP results) and the SLD landmark output; they are written by the last `-Fresh` run.
The reported results come from a `-Fresh` run, and a second `-Fresh` run reproduced all
predictions byte for byte.

Across computers the GCZ geometry is not guaranteed to be identical. The previous version of this
package was checkpointed on another computer with the same code and the same package versions;
there, the ICP alignment of atlas S03 L3 to S02 L3 reached a residual of 1.33 mm, whereas on the
computer used for the reported results it converged to a worse fit and atlas S01 L1 (1.51 mm) was
selected instead. Such floating-point differences change some atlas choices and therefore the
training samples and the predictions (in that earlier run: GCZ 3.208 mm, SLD 2.224 mm, and GCZ-RL
without the side constraint 3.336 mm). The checkpoints in `checkpoints/` are those of the
reported run, so `RUN_ALL.ps1 -Part 2` without `-Fresh` reproduces the reported numbers on any
computer.

## 2. Tables and figures → code → output file

The items below are the output files of this part. Their place in the journal manuscript
(Tables 4 and 5, the leakage test of Section 4.3, Supplementary Table S9) is listed in `docs/PAPER_MAP.md`.

| Item | What it shows | Produced by | Output file (`outputs/`) |
|---|---|---|---|
| Table P2-1 | Folds: training/test subjects and vertebrae, training samples per ligament | `run_part2.py`, `run_fold()` | `table_P2_1_folds.csv` |
| Table P2-2 | Overall MAE, Pass@5 and coverage metrics of GCZ, GCZ-RL, SLD | `run_part2.py`, `main()` | `table_P2_2_overall.csv` |
| Table P2-3 | MAE and Pass@5 per ligament | `run_part2.py`, `main()` | `table_P2_3_per_ligament.csv` |
| Table P2-4 | MAE and Pass@5 per held-out subject | `run_part2.py`, `main()` | `table_P2_4_per_subject.csv` |
| Table P2-5 | Data-leakage test | `run_leakage_test.py` | `table_P2_5_leakage_test.json` |
| Table P2-6 | Extra Trees configuration and α selected in each fold | `run_part2.py`, `run_fold()` | `table_P2_6_selected_models.csv` |
| Figure P2-1 | Per-ligament MAE of the three methods | `run_part2.py`, `figure()` | `figure_P2_1_per_ligament.png/.pdf` |

Also in `outputs/`: predicted coordinates (`predictions_gcz.csv`, `predictions_gczrl.csv`,
`predictions_sld.csv`), per-vertebra metrics, the complete inner selection grid
(`inner_selection_grid.csv`), the trained models of every fold (`models/`), fold details
(`fold_details.json`) and the run log.

## 3. Cohort and split

* All annotated lumbar vertebrae of the dataset: **S01 L1–L5 (5), S02 L1–L4 (4), S03 L1–L5 (5) =
  14 vertebrae from 3 subjects**. S02 has four lumbar vertebrae because its scan ends at L4.
  S04 and S05 are annotated cervical subjects (C3–C7) with no lumbar vertebrae.
* Three folds; each holds out one complete subject:

  | Held-out subject | Training subjects | Training vertebrae | Test vertebrae | Training samples |
  |---|---|---|---|---|
  | S01 | S02 + S03 | 4 + 5 = 9 | 5 | 741 |
  | S02 | S01 + S03 | 5 + 5 = 10 | 4 | 836 |
  | S03 | S01 + S02 | 5 + 4 = 9 | 5 | 765 |

## 4. How GCZ-RL (Extra Trees) is trained with 14 vertebrae

The model never sees all 14 vertebrae. In each fold it learns from the 9 or 10 vertebrae of
the two training subjects and is applied to the third subject.

1. **Zones and atlas** come from the training subjects' annotations and meshes.
2. **Training-side GCZ.** Each training vertebra is predicted by GCZ using the zones and atlas of
   the *other* training subject, exactly as a new patient would be predicted.
3. **Samples.** The training unit is a landmark point, not a vertebra. GCZ gives about 10 points
   per ligament; each point is paired with an expert landmark and its residual (expert − GCZ)
   becomes a target; bilateral ligaments are mirrored. This gives 741–836 samples × 36 features
   per fold. Per ligament (fold S01): ALL 90, PLL 90, LF 180, ISL 164, CL 172, ITL 36, SSL 9.
4. **Regressor.** Only Extra Trees is used: four configurations (ET-1 … ET-4, the same four that
   form the dissertation's final model; definitions in `code/ligament/gczrl.py`). For every
   ligament, the configuration and the blend weight α ∈ {0, 0.25, 0.40, 0.55, 0.70, 0.85, 1.00,
   1.10, 1.20} are chosen by an **inner leave-one-subject-out** on the two training subjects
   (train on one, evaluate on the other, both directions). α = 0 means that the inner validation
   found no benefit from the correction, so the GCZ point is kept.
5. **Final fit.** The selected Extra Trees are trained on all samples of the two training subjects
   and applied once to the held-out subject.
6. **Side constraint for lateral ligaments (ITL, CL).** After correction, each point is projected
   back onto the ligament zone. For ITL and CL the projection only uses zone vertices on the same
   side of the vertebral midline (x-centroid of the mesh, the midline also used for mirroring) as
   the GCZ point being corrected; if the zone has no vertex on that side, the correction is not
   applied and the GCZ point is kept (`SIDE_LOCKED` and `project_to_zone` in
   `code/ligament/gczrl.py`). The rule uses only the mesh and the GCZ point, never annotations, and
   is applied in the inner selection as well as in the final prediction. It was added after a run
   without it placed both ITL points of S01 L5 and S02 L4 on the same transverse process
   (ITL MAE about 50 mm in those two vertebrae). It is switched off in Part 1, which reproduces
   the dissertation unchanged.

In Part 1 (dissertation) the same method was trained once on the 20 calibration vertebrae
(1683 samples), with the configuration and α of each ligament fixed in
`FINAL_CONFIGURATION`.

## 5. Data-leakage controls

* The held-out subject's annotations are not reachable during training: the label function
  returns nothing for test vertebrae (`labels` in `run_fold()`); they are read only in the final
  evaluation.
* Zones, atlas, benchmark prior, correction radii, regressor configuration and α are computed
  from the two training subjects. Inside the inner selection, zones and radii are recomputed
  without the inner held-out subject, and mirrored rows are not evaluated.
* **Leakage test (Table P2-5).** The S02 fold is recomputed from scratch twice, once with the true
  S02 annotations and once with every S02 landmark corrupted by 25 mm Gaussian noise. All 424
  predicted S02 points of the two runs are identical (maximum difference 0 mm). Both runs are
  recomputed, so the test does not depend on the checkpoints; it exits with an error if any
  point differs.
* `code/verify_all.py` checks that every saved model of every fold is an
  `ExtraTreesRegressor` and that no training vertebra belongs to the held-out subject.

## 6. Results

| Method | MAE (mm) | Pass@5 (%) | Reference coverage (%) | Missing outputs |
|---|---|---|---|---|
| GCZ | 2.879 | 86.0 | 58.4 | 0 |
| GCZ-RL (Extra Trees) | 2.150 | 92.7 | 58.4 | 0 |
| SLD | 2.224 | 90.2 | 67.0 | 0 |

Per ligament (MAE mm / Pass@5 %, Pass@5 pooled over the matched points of each class):

| | ALL | PLL | LF | ISL | SSL | ITL | CL |
|---|---|---|---|---|---|---|---|
| GCZ | 2.61 / 92.1 | 1.86 / 100.0 | 1.92 / 97.9 | 2.29 / 87.1 | 2.92 / 78.6 | 6.26 / 53.6 | 2.30 / 92.9 |
| GCZ-RL | 2.38 / 97.1 | 1.86 / 100.0 | 1.91 / 97.9 | 2.23 / 87.9 | 2.69 / 78.6 | 1.80 / 92.9 | 2.18 / 94.3 |
| SLD | 1.66 / 99.5 | 2.45 / 88.8 | 5.00 / 58.3 | 2.04 / 93.4 | 1.97 / 92.9 | 0.87 / 100.0 | 1.59 / 98.9 |

Per held-out subject (MAE mm): GCZ 3.78 / 2.21 / 2.51, GCZ-RL 2.32 / 1.99 / 2.12,
SLD 1.77 / 2.26 / 2.65 for S01 / S02 / S03.

The committee report (`committee_results_report_neutral.pdf`, built by `report/build_report.py`)
uses the same predictions; its Pass@5 table averages per vertebra–ligament case instead of pooling
points, which changes a few per-ligament values by at most 1.2 percentage points.

Interpretation:

* GCZ-RL improves on GCZ for every held-out subject (MAE −38.8 %, −10.1 %, −15.7 % for S01, S02,
  S03), for 13 of 14 vertebrae, and for 6 of 7 ligament classes (PLL unchanged); Pass@5 rises
  from 86.0 to 92.7 %.
* GCZ-RL and SLD are comparable overall (MAE 2.150 vs 2.224 mm, Pass@5 92.7 vs 90.2 %). SLD has
  the lower MAE for 5 of 7 classes and for 8 of 14 vertebrae; GCZ-RL is lower for LF and PLL,
  where SLD is weakest. The overall difference is small and is not a demonstrated advantage.
* ITL: with the side constraint, GCZ-RL's ITL error in S01 L5 and S02 L4 is 4.96 and 5.53 mm (the
  ITL zone has no vertex on the left side there, so the left point keeps its GCZ position); on the
  other 12 vertebrae it is 0.6–2.3 mm.
* In 7 of the 21 ligament–fold choices the inner validation selected α = 0, so the correction is
  only used where the two training subjects support it.
* With three subjects, these are descriptive results; no confidence intervals or significance
  tests are given.

## 7. SLD in this part

SLD is run with the official SpineToolkit/SpineLib code (`code/third_party/SpineLib`, unmodified)
and its official preprocessing, through `code/ligament/sld_official.py`:

* the meshes are stored in LPS coordinates and are converted to the RAS frame that SpineLib uses
  (as 3D Slicer does when loading them), and the landmarks are converted back for evaluation;
* 3D Slicer objects are replaced by a small stub (`code/third_party/slicer_stub`); the VMTK
  centerline step, which is not available without Slicer, is approximated by cross-section
  centroids along a surface path.
SLD has no trainable parameters and uses no annotations. These SLD values differ from the SLD
column of dissertation Table 4.4, which came from an earlier standalone run without the
coordinate-frame conversion.

## 8. Limitations

* The experimental data are what was available: three annotated subjects (14 lumbar vertebrae)
  from one public CT dataset, annotated once. Training folds contain only two subjects.
* Ligaments with one or two points per vertebra (SSL, ITL) have very few training samples
  (9–40 per fold), which limits what the residual regressor can learn for them.
* Matched-point metrics do not penalise unmatched reference points; coverage metrics are given
  in Table P2-2.
* The ITL/CL side constraint (Section 4, step 6) was introduced after the held-out results of a
  run without it had been seen. It is label-free and anatomically motivated, but it is a post-hoc
  correction and is reported as such.
* SLD's VMTK centerline step is approximated (Section 7); the official Slicer module may perform
  differently.
* GCZ's ICP atlas alignment is sensitive to floating-point differences between computers
  (Section 1), so a `-Fresh` run on another computer can give different GCZ and GCZ-RL values.
* No inter-observer data, no pathology, and no evaluation of the effect on finite element results.
