# GCZ-RL: Geometry-Calibrated Residual Learning for Spinal Ligament Landmark Localization on Vertebral Surface Meshes

Code, trained models, predictions, and reproducibility scripts for the manuscript

> M. Ahmadi, M. Lin, Y. Tang, F. D. Vrionis. *GCZ-RL: Geometry-Calibrated Residual Learning for
> Spinal Ligament Landmark Localization on Vertebral Surface Meshes.*
> 2026.

GCZ-RL localizes sparse ligament attachment landmarks (ALL, PLL, LF, ISL, SSL, ITL, CL) on
vertebral surface meshes. It works in three steps:

1. **Calibrated regions.** It learns ligament-specific admissible regions from expert
   annotations in normalized vertebral coordinates (geometry-calibrated zoning, GCZ).
2. **Atlas transfer.** It transfers an annotated atlas by rigid ICP and snaps the transferred
   landmarks into these regions.
3. **Bounded correction.** It corrects each landmark with a bounded Extra Trees residual and
   projects it back onto its region.

Every final prediction is a vertex of the target surface. Corrected landmarks are projected into
the admissible region of their class (740 of 742 points in the primary evaluation); the other two
are ITL points that keep their GCZ position under the side-preserving rule, because the region has
no vertex on their side.

![Graphical abstract](paper/outputs/graphical_abstract.png)

## Main results

Primary evaluation: patient-level leave-one-subject-out on all 14 lumbar vertebrae of three
subjects. Calibration, atlas, residual targets, correction bounds and model selection are all
nested inside the training subjects.

| Method | MAE (mm) | Pass@5 (%) | Closest point P→R (mm) |
|---|---|---|---|
| GCZ (geometry only) | 2.879 | 86.0 | 2.873 |
| SLD (released code, same evaluator) | 2.224 | 90.4 | 1.787 |
| **GCZ-RL** | **2.150** | **92.7** | 2.140 |

- **GCZ-RL vs GCZ.** Difference −0.729 mm, 95% vertebra-level bootstrap CI −1.078 to −0.381;
  Wilcoxon p < 0.001. Every held-out subject improved.
- **GCZ-RL vs SLD.** Difference −0.074 mm, CI −0.386 to 0.228; p = 0.35, not significant.

MAE is the mean Euclidean distance after one-to-one Hungarian matching. Pass@5 is averaged
over vertebra–class cases.

A supplementary vertebra-wise development analysis (57 vertebrae, cervical to lumbar) is in
`part1_dissertation_57/`.

## Repository layout

```
code/ligament/            method and evaluation (gcz.py, gczrl.py, metrics.py, data_io.py, sld_official.py)
code/third_party/         3D Slicer stub for running SLD outside Slicer; SpineLib is fetched, not stored
code/verify_all.py        independent re-computation of the reported metrics (9 checks)
dataset/                  split, subject codes, checksums, source manifest; raw data are fetched
part1_dissertation_57/    development analysis (20 calibration / 37 held-out vertebrae)
part2_patient_loso_14/    primary analysis: patient-level LOSO, leakage test
paper/                    scripts for every number, table and figure added for the journal paper
scripts/                  fetch_data.py, fetch_spinelib.py
docs/                     PAPER_MAP.md (paper item -> script -> output), CODE_PAPER_AUDIT.md
tests/                    unit tests of the geometric constraints and the metrics
manuscript/               LaTeX source and PDFs of the paper, cover letter, Statement of Significance
```

## Installation

Python 3.11 is required, and the versions in `requirements.txt` are pinned.

```bash
python -m venv .venv
.venv/bin/pip install -r requirements.txt      # Windows: .venv\Scripts\pip install -r requirements.txt
```

## Data

The vertebral meshes and expert ligament annotations were published by the SLD authors in the
[SpineToolkit](https://github.com/VisSim-UniKO/SpineToolkit) repository, and are based on VerSe CT
data. That repository has no licence file, so this project **does not redistribute** them.
Instead, one command downloads them from a pinned commit, rebuilds the exact input layout, and
verifies all 1678 files against the checksums of the files used for the paper:

```bash
python scripts/fetch_data.py                        # -> dataset/raw (about 200 MB)
python scripts/fetch_data.py --out D:/gczrl_raw     # elsewhere; then set GCZRL_RAW_ROOT=D:/gczrl_raw
```

On Windows, use `--out` with a short path if the repository sits in a deep folder (the 260-character
path limit). See [dataset/README.md](dataset/README.md) for the subject mapping and the provenance of every file.

## Reproduce the paper

One command fetches the data if needed, runs every analysis, regenerates every figure and table
of the paper, checks the manuscript against the regenerated outputs, and runs the unit tests:

```powershell
powershell -ExecutionPolicy Bypass -File RUN_ALL.ps1 -Pdf    # Windows; -Pdf also compiles the manuscript (needs latexmk, e.g. MiKTeX)
```

```bash
PDF=1 bash run_all.sh                                       # Linux and macOS
```

The run takes about an hour. Add `-Slow` (`SLOW=1`) to also recompute Table 3 and
Supplementary Table S6 (about 30 min), the original notebook (Supplementary Figs. S1–S2,
Tables S7–S8), and the ten ICP seed runs of Fig. 10 and Supplementary Table S16 (2–3 hours);
without it, these are read from their saved outputs. The original notebook also writes about
1.5 GB of interactive HTML figures, which are not part of the repository; keep at least 3 GB free.

Then compare every regenerated output with the committed version:

```bash
python scripts/compare_outputs.py      # Windows: .venv\Scripts\python.exe scripts\compare_outputs.py
```

It classifies every changed file as same data (CSV and JSON values within 1e-9), rendering noise
(images, and PDF pages rasterized at 100 dpi), provenance (logs, run times, SVG dates, and the
absolute file paths that the original notebook prints), or CHANGED, and exits with an error if
anything is CHANGED.

The individual steps, in the order of the run scripts:

```bash
python part1_dissertation_57/run_part1.py          # development analysis (Table 9)
python part2_patient_loso_14/run_part2.py          # primary LOSO (Tables 4-5), minutes with stored ICP results
python part2_patient_loso_14/run_leakage_test.py   # data-leakage test
python code/verify_all.py                          # independent checks -> VERIFICATION_REPORT.json
python paper/make_paper_numbers.py                 # statistics, Table 6, Supplementary S2-S4; checks every quoted value
python paper/fig6_qualitative_loso.py              # Supplementary Fig. S7 (the previous qualitative figure)
python paper/check_fig6.py                         # Fig. 7 and S7 points = saved predictions; caption numbers
python paper/graphical_abstract.py                 # graphical abstract
python paper/revision_analyses.py                  # component analyses (Section 4.4, Supplementary S8-S11)
python paper/figures_main.py                       # Figs. 6, 9, 10, S3-S5 -> paper/outputs/figures (audit CSV)
python paper/region_verification.py                # Supplementary Table S17
python paper/margin_verification.py                # Supplementary S16, Fig. S6
python paper/figures_anatomy.py                    # Figs. 7, 8, S8
python paper/fig11_development.py                  # Fig. 11
python paper/itl_trace.py                          # Supplementary S19, Table S18
python paper/itl_double_check.py                   # independent checks of the ITL result (Supplementary S19)
python paper/icp_seed_sensitivity.py summary       # Table S16 from the saved runs ("run 0 1 ... 9" recomputes them)
python paper/revision_tables.py                    # table rows of Tables 6-8 and S10-S14 -> manuscript/tables
python paper/reproducibility_check.py              # 68 reported claims vs regenerated values (Table S15)
python paper/check_text_numbers.py                 # numbers of the revised text vs regenerated values
python paper/check_tables.py                       # every cell and bold marking of the typed tables
python paper/build_manuscript_assets.py --pdf      # figures and tables into manuscript/; compiles the PDFs
python -m unittest discover -s tests               # geometric constraints and metrics
python paper/table3_regressor_comparison.py        # Table 3 and Supplementary Table S6 (slow step, about 25 min)
```

`--fresh` recomputes the stored ICP geometry and, for Part 2, the SLD output. SLD needs the
SpineLib files, which you fetch with `python scripts/fetch_spinelib.py`.

> **Reproducibility note.** With the stored ICP alignments (the default), every number in the
> paper is regenerated exactly from the public data. This was verified from a clean download.
> Recomputing the alignments with `--fresh` is deterministic within one environment, but it can
> select a different atlas vertebra for individual targets in a different environment.
>
> In a second environment, 2 of the 14 held-out vertebrae changed atlas. The result was
> GCZ 3.208 mm and GCZ-RL 2.418 mm, against 2.879 and 2.150 mm in the paper; SLD was
> unchanged. See Supplementary Section S5 and [docs/CODE_PAPER_AUDIT.md](docs/CODE_PAPER_AUDIT.md).
>
> `paper/icp_seed_sensitivity.py` quantifies this: it re-runs the complete nested evaluation with
> all ICP alignments recomputed for ten seed offsets, with nothing else changed. GCZ-RL improved on
> GCZ in every run (22.8–33.4 %, mean 26.9 %) and for every held-out subject in every run; its
> absolute error ranged from 2.10 to 2.92 mm (mean 2.36 mm). Offset 0 reproduces the
> second-environment result exactly (Supplementary Section S13).

[docs/PAPER_MAP.md](docs/PAPER_MAP.md) lists, for every table and figure of the paper, the
script that produces it and its output file.

## Component analyses

`paper/revision_analyses.py` examines each part of the method on the primary evaluation, using
the stored ICP geometry and the fold models written by `run_part2.py`. It selects nothing, and
its complete configuration reproduces the saved GCZ-RL predictions exactly.

| Configuration | MAE (mm) |
|---|---|
| GCZ (geometry only) | 2.879 |
| GCZ + projection onto the calibrated regions (strongest geometry-only configuration) | 2.310 |
| GCZ + bounded residual, no projection | 2.698 |
| GCZ + bounded residual + projection | 2.132 |
| **GCZ-RL (complete, nested blend weights)** | **2.150** |

- **Strongest geometry-only configuration.** GCZ-RL improves on it for each of the three
  held-out subjects; over all cases the difference is −0.160 mm (95% vertebra-level bootstrap CI
  −0.234 to −0.085; Wilcoxon p < 0.001 over the 98 vertebra–class cases of this cohort).
- **Candidate regions.** The regions contain 90.0% of the held-out reference landmarks while
  covering a median of 0.1–26.9% of the vertebral surface area.
- **Learned residuals.** They point in the direction of the true residual for 67.1% of
  held-out points.
- **Atlas choice.** The ICP residual of the selected atlas is not associated with landmark
  error (Spearman ρ = 0.24).

Outputs are in `paper/outputs/revision/`. `paper/reproducibility_check.py` verifies 68 main
claims of the paper against the regenerated outputs.

## Notes on the evaluation

- **Matched-point metrics.** The reference annotations are dense point sets, while GCZ and
  GCZ-RL output at most ten landmarks per class. MAE and Pass@5 therefore score matched points
  only. Point coverage and matching-independent metrics are reported alongside them
  (Supplementary Tables S2–S3).
- **SLD.** SLD is run from its released SpineLib code outside 3D Slicer. The VMTK centerline is
  replaced by a documented cross-section approximation
  (`code/third_party/slicer_stub/ExtractCenterline.py`). The SLD numbers refer to this
  implementation.
- **Post hoc constraint.** The left/right side constraint for ITL and CL in the primary analysis
  was introduced after an earlier run had been inspected; this is disclosed in the paper.
- **Audit.** [docs/CODE_PAPER_AUDIT.md](docs/CODE_PAPER_AUDIT.md) records the code-versus-manuscript
  audit done before submission, including one corrected bug in the development regressor
  comparison.

## Citation

See [CITATION.cff](CITATION.cff). Please also cite SLD (Blomenkamp et al., ISBI 2026) and VerSe
(Sekuboyina et al., Med. Image Anal. 2021) when you use the data.

## Licence

The code in this repository is released under the MIT licence ([LICENSE](LICENSE)). The licence
does not cover the source data or SpineLib, which remain under the terms of their owners.
