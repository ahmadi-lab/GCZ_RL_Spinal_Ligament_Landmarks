# Code ↔ manuscript audit (October 2026)

This audit was run before the *Medical Image Analysis* submission. Every number, equation, and
method statement of the manuscript was compared with the code and the saved outputs of this
repository. The table numbers below refer to the revised manuscript.

## 1. Verified without change

| Item | Result |
|---|---|
| Primary LOSO (Tables 4–5) | All per-ligament and per-subject MAE and Pass@5 values were regenerated exactly. This covers the macro values 2.879 / 2.150 / 2.224 mm and 86.0 / 92.7 / 90.4 %, the reductions of 25.3 % and 38.8 / 10.1 / 15.7 %, and the 741 / 836 / 765 training samples. |
| Text claims, primary | GCZ-RL beat GCZ for 13 of the 14 vertebrae and for 6 of the 7 classes (PLL unchanged). The two side-empty ITL cases have errors of 4.96 and 5.53 mm; the other twelve lie between 0.6 and 2.3 mm. α = 0 was selected in 7 of 21 class–fold choices. The leakage test left all 424 points unchanged. |
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
| `run_part2.py` (stored ICP results) | see `VERIFICATION_REPORT.json` |
| `run_part1.py` (stored ICP results) | see `VERIFICATION_REPORT.json` |
| `code/verify_all.py` | 9 / 9 checks |
| `paper/make_paper_numbers.py` | all manuscript values reproduced |

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
