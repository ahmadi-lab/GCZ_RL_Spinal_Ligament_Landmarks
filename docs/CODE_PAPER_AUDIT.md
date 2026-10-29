# Code ↔ manuscript audit (October 2026; last updated 2026-10-29)

This audit was run before the journal submission. Every number, equation, and
method statement of the manuscript was compared with the code and the saved outputs of this
repository. The table numbers below refer to the revised manuscript.

## 1. Verified without change

| Item | Result |
|---|---|
| Primary LOSO (Tables 4–5) | All per-ligament and per-subject MAE and Pass@5 values were regenerated exactly. This covers the macro values 2.879 / 2.150 / 2.224 mm and 86.0 / 92.7 / 90.4 %, the reductions of 25.3 % and 38.8 / 10.1 / 15.7 %, and the 741 / 836 / 765 training samples. |
| Text claims, primary | GCZ-RL had the lower error for 13 of the 14 vertebrae (S03 L1: 0.076 mm higher) and for 6 of the 7 classes (PLL unchanged). The two side-empty ITL cases have errors of 4.96 and 5.53 mm; the other twelve lie between 0.6 and 2.3 mm. α = 0 was selected in 7 of 21 class–fold choices. The leakage test left all 424 points unchanged. |
| Supplementary S2–S4 | Point accounting, matching-independent metrics, and region feasibility all match. Three region-size medians were truncated in the paper (for example 5254 instead of 5254.5); they are now printed exactly. |
| Development (Table 6) | All 24 values match. The regenerated GCZ and GCZ-RL predictions are identical to the archived ones (259 / 259 point sets). |
| Features (Table S1) | All 36 features match `point_features` and `_zone_features`. |
| Matching diagnostic (Tables S7–S8) | All 44 values match on a re-run of the original notebook (`NB/outputs/tables/check_against_paper.csv`). |
| Data | Every annotation used is byte-identical to the public SpineToolkit release (see §3). |

## 2. Discrepancies found and how they were resolved

| # | Finding | Resolution |
|---|---|---|
| 1 | **Correction equations.** The code projects twice: once after the clipped correction, and again after the α-blend (`corrected_points` → `blend`). The paper described one projection. Because α reaches 1.2, the cap is α·τ, not τ. | The manuscript equations (Eqs. 9–10) were rewritten to match the code. |
| 2 | **Mirroring.** Mirroring augments LF, ISL, ITL, and CL about the sagittal plane through the vertex centroid. The paper said "bilateral classes" mirrored about "the mid-sagittal plane", but ISL is a midline structure. | The manuscript now names the classes and the plane. |
| 3 | **GCZ selection.** The paper did not mention the 15 mm snap radius, the side and quadrant snapping, or the farthest-point sampling to ten landmarks. It also did not say that residual targets use Hungarian pairing. | All added to Sections 3.6–3.7 and Table S5. |
| 4 | **Bug in the calibration regressor comparison (Table 3).** `held_out_records(original_rows_only=False)` validated the first *n* rows of each bilateral class. These rows interleave original and mirrored samples, so half of the predicted corrections were applied to the wrong GCZ points. Zones, class prior, and correction bound also included the held-out vertebra. | Re-run with `original_rows_only=True` and per-fold quantities (`paper/table3_regressor_comparison.py`). Extra Trees now has the lowest error for **all 7 classes** (mean 2.037 mm; random forest 2.307, gradient boosting 2.337). Before the fix, random forest was lowest for LF. The old values are given in Table S6. The final model choice is unchanged. |
| 5 | **Development baseline (Table 6).** The "GCZ" column and the GCZ points that GCZ-RL corrects come from two different runs: STL against VTK surfaces, different seeds, and a different last fallback. | The manuscript discloses this and reports the corrected input (3.704 mm, Pass@5 80.8 %). A paired test against that input was added (Wilcoxon p < 0.001). |
| 6 | **Matching diagnostic.** The notebook behind the old Tables 3–4 sets the number of GCZ landmarks per class to the number of reference points of the evaluated vertebra, which uses annotation cardinality. | The analysis was moved to Supplementary S6 with an explicit disclosure. It is not used for any method comparison. |
| 7 | **Alignment table.** The table did not say which vertebrae it covers (the 9 development lumbar vertebrae), and there was no alignment check for the primary analysis. | Both are now reported: primary atlas residual 1.56 mm (range 1.14–2.16); the same lumbar level was chosen in 7 of 14 cases. |
| 8 | **Pass@5 conventions.** The case-averaged value is 90.4 % for SLD and the class-pooled value is 90.2 %. The committee README reported 90.2. | The manuscript states both. For GCZ and GCZ-RL the two conventions coincide. |
| 9 | **SLD's published error.** SLD reports 0.7 mm with a closest-point error; this evaluation gives 2.224 mm. | Closest-point metrics in both directions were added. Under SLD's own definition SLD leads (1.787 vs 2.140 mm), and the Discussion explains the gap. |
| 10 | **No statistical tests.** | Wilcoxon signed-rank tests and a vertebra-level cluster bootstrap were added. GCZ-RL vs GCZ: −0.729 mm [−1.078, −0.381]. GCZ-RL vs SLD: −0.074 mm [−0.386, 0.228], p = 0.35. |
| 11 | **References.** SLD is now published (ISBI 2026, DOI 10.1109/ISBI61048.2026.11515825). The Biswas preprint is 2025, not 2026. *World Neurosurg.* is volume 201. Some DOIs were missing. | All 27 DOIs were resolved against Crossref, and author names come from the Crossref records. |
| 12 | **Highlights.** "Improves … for every class" is false, because PLL is unchanged. "Official SLD" overstates the comparison, because the centerline is substituted. | The highlights were rewritten. |
| 13 | **CRediT.** "Clinical guidance" is not a CRediT role. | Flagged for the authors. |
| 14 | **ICP geometry not reproducible across environments.** Recomputing the primary-evaluation alignments from scratch here, with the same code, pinned library versions, and checksum-verified data, changed the atlas of 2 of 14 held-out vertebrae (S01 L3, S02 L3) and the inner model selection of the S01 and S02 folds. The result was GCZ 3.208 mm / 83.1 % and GCZ-RL 2.418 mm / 89.8 %, against 2.879 / 2.150 mm in the paper; SLD and the S03 fold were unchanged. The cause is that atlas selection runs ICP from raw scanner coordinates. Many candidates converge to wrong minima (4–5 mm residuals), and for S02 L3 the lower-residual atlas (1.33 vs 1.51 mm) gave a 4.0 mm worse GCZ error. Forcing single-threaded BLAS does not change this, and the run is deterministic within this environment. | **Disclosed** in the main-text limitations and Supplementary S5. The stored checkpoints regenerate the paper exactly. **Recommended:** make atlas selection robust (centroid or principal-axis pre-alignment, multi-start ICP, orientation check), re-run all analyses, and report variability over seeds. |

## 3. Data provenance

- **Annotations.** All 1678 input files can be rebuilt from SpineToolkit commit `62f39cf`
  (`scripts/fetch_data.py`, `dataset/source_manifest.csv`). The 759 VTK-folder files and the S04
  files are byte-identical to the official ones. The 712 annotation copies in the `*_STL` folders
  are identical to the official annotations.
- **STL meshes.** The STL files used for the development GCZ column are not part of the official
  release. They are rebuilt from the official VTK meshes with identical vertex arrays. In
  `Verse581_T1`, one degenerate triangle has reversed corner order; the fetch script reproduces it.
- **Frame mismatch.** The S04 STL files declare `SPACE=RAS` while their annotations declare LPS.
  S04–S05 are never scored. Check this before ever scoring them.
- **No licence.** SpineToolkit has no licence file. The data and SpineLib are therefore fetched,
  not redistributed. The SpineLib copy used in the analyses was an unmodified subset of 21 files,
  verified by git blob SHA-1.

## 4. Reproduction from public sources

On 2026-10-21 the data were fetched from GitHub into an empty folder (1678 / 1678 files
verified). The analyses were then re-run from this repository:

| Run | Result |
|---|---|
| `scripts/fetch_data.py` | 1678 / 1678 files verified |
| `run_part2.py` (stored ICP results) | Every table, prediction file, and the inner selection grid is **byte-identical** to the committed outputs. |
| `run_part1.py` (stored ICP results) | Table 6 and all 259 prediction sets are identical. The regressor grid differed only in the third decimal of gradient boosting, from multi-threaded fitting; histogram gradient boosting is now fitted with one thread, and the grid is reproduced exactly (2026-10-28). |
| `run_leakage_test.py` | PASS: 0 mm difference. It compares 412 points here against 424 in the original environment, because it recomputes the ICP geometry (finding 14). |
| `code/verify_all.py` | 9 / 9 checks passed |
| `paper/make_paper_numbers.py` | 46 / 46 manuscript values reproduced |
| `run_part2.py`, ICP recomputed | GCZ 3.208, GCZ-RL 2.418, SLD 2.224 mm (finding 14) |

## 5. Open items for the authors

1. **SLD centerline.** SLD is compared with a substituted centerline. A run inside 3D Slicer with
   VMTK, ideally with the SLD authors, would make the comparison definitive.
2. **External test set.** No external test set is used. The SpineToolkit release also contains
   annotated artificial Sawbone vertebrae (C1–S1), a possible external, if synthetic, test.
3. **Midplane values.** The midplane values in Fig. 2D come from a figure of the original
   notebook (cells 21–24) and are not regenerated as numbers.
4. **Paths in the original notebook.** The original notebook keeps the local absolute paths of
   the machine it was run on, as a record of provenance. Decide whether to keep them in the
   public release.

## 6. Verification for the revision (2026-10-27)

Before the manuscript was revised, the data were fetched again (1678 / 1678 files verified) and
every analysis was re-run from commit `17d5661` with the pinned environment.

| Run | Result |
|---|---|
| `paper/make_paper_numbers.py` | 46 / 46 manuscript values reproduced |
| `run_part2.py` (stored ICP results) | Predictions, metrics, and selections byte-identical to the committed outputs |
| `run_part1.py` (stored ICP results) | Table 4.4 exact; 259 / 259 point sets identical; final-model checksum unchanged |
| `paper/table3_regressor_comparison.py` | All values identical |
| `code/verify_all.py` | 9 / 9 checks passed |
| `run_leakage_test.py` | PASS (412 points, as in finding 14) |
| `paper/reproducibility_check.py` | 68 / 68 claims verified; the S02 reduction is 10.0 %, not 10.1 % (it had been computed from rounded means) |

The component analyses of the revision (`paper/revision_analyses.py`) reproduce the saved
GCZ-RL predictions from the stored geometry and fold models with a maximum coordinate difference
of 6e-14 mm. Files that differ between runs only in provenance (time stamps, `"gcz_geometry":
"checkpoint"`) were kept at their committed versions.

### ICP seed sensitivity (2026-10-27)

`paper/icp_seed_sensitivity.py` re-ran the complete nested evaluation ten times with every ICP
alignment recomputed and the ICP sampling seeds of `run_part2.gcz_predict` shifted by
10000 × offset. Folds, calibration, fallback and farthest-point-sampling seeds, Extra Trees
seeds, selection grids, and evaluation were unchanged; the stored checkpoints were not touched.
Offset 0 reproduced the second-environment result of finding 14 exactly (3.208 / 2.418 mm).
Over the ten offsets, 2–4 of 14 atlases differed from the stored alignments; GCZ-RL improved on
GCZ in every run (22.8–33.4 %) and for every held-out subject in every run; GCZ-RL ranged from
2.103 to 2.918 mm (mean 2.357, SD 0.239) and GCZ from 2.844 to 3.888 mm. GCZ-RL was below SLD
(2.224 mm) in 3 of 10 runs.

## 7. Final verification from a fresh clone (2026-10-29)

Commit `3906add` was cloned into an empty folder and run on Windows 11 with the pinned environment
(Python 3.11.15, `requirements.txt`) and the public data (`scripts/fetch_data.py`, 1678 / 1678
files verified): first `RUN_ALL.ps1 -Pdf` (every analysis, figure, table, manuscript check, the
compiled PDFs, and the unit tests; 60 min), then `RUN_ALL.ps1 -Part 1 -Slow` (Part 1 again, the
original notebook, Table 3 and Supplementary Table S6, and the ten ICP seed runs, all recomputed;
4 h, of which 3 h 6 min for the ICP runs). Every step finished without error.

| Check | Result |
|---|---|
| `code/verify_all.py` | 9 / 9 checks passed |
| `paper/make_paper_numbers.py` | 46 / 46 manuscript values reproduced |
| `run_part1.py` | Table 4.4 exact; 259 / 259 GCZ-RL point sets identical |
| `run_leakage_test.py` | S02 predictions unchanged with corrupted S02 annotations (0 mm, 412 points with ICP recomputed) |
| `paper/check_fig6.py` | Drawn points = saved predictions; caption errors 3.263 / 1.972 / 2.771 mm |
| `paper/reproducibility_check.py` | 68 / 68 reported values |
| `paper/check_text_numbers.py` | 55 / 55 passages of the revision |
| `paper/check_tables.py` | 512 / 512 table cells, 35 / 35 bold markings |
| `paper/itl_double_check.py` | 12 / 12 checks passed |
| Unit tests | 10 / 10 |
| Original notebook | 37 code cells without error, 32 with output identical to the saved notebook (the other 5 differ in file paths, library versions, or a figure saved as HTML); 44 / 44 published values |
| `paper/table3_regressor_comparison.py` | All four output tables byte-identical |
| `paper/icp_seed_sensitivity.py run 0 … 9`, `summary` | All ten runs identical to the committed results except the recorded run time; summary byte-identical |

`scripts/compare_outputs.py` then compared every output with its committed version. The two runs
rewrote 70 committed files, and none changed: 35 PDFs (figures, graphical abstract, and the clean
and highlighted manuscript and the supplement) were pixel-identical on every page; three CSV files of
the ablation and interpretability analyses differed by at most 8e-15 (last floating-point digits of tree averages computed
in parallel); the other 32 differed only in run provenance (run times of the ten ICP runs and of the
notebook cells, time stamps of the run logs, SVG creation dates, and the absolute file paths that
the original notebook prints). Every other regenerated file, every PNG included, was byte-identical.

Changes after `3906add`, none of which changes an analysis or an output: `scripts/compare_outputs.py`
classifies absolute file paths and notebook cell run times as provenance, matches provenance keys as
whole words (so that a key such as "Candidate-region relation" is compared), and lists new files;
`RUN_ALL.ps1` no longer writes a scikit-learn warning that is repeated for every nested parallel
Extra Trees fit (over 1 GB of log with `-Slow`) and logs only its count; the run times stated in the
run scripts and the README were corrected.
