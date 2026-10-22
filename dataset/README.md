# Data

## Sources

| Subjects | Vertebrae | Meshes | Expert ligament annotations |
|---|---|---|---|
| S01–S03 (VerSe `verse005`, `verse561`, `verse581`) | 47 (C3–L5) | SpineToolkit `VerSe_Ligaments_GT` (VTK), derived from VerSe CT | SpineToolkit `VerSe_Ligaments_GT` |
| S04–S05 (`Bc11701`, `Bc11702`) | 10 (C3–C7) | SpineToolkit `Cervical_Ligaments_GT` | SpineToolkit `Cervical_Ligaments_GT` |

- **Source repository.** Both annotation datasets were published by Blomenkamp et al. with SLD
  (ISBI 2026) at <https://github.com/VisSim-UniKO/SpineToolkit>, folder
  `LigamentLandmarkDetection/Datasets`, commit `62f39cf`.
- **VerSe.** The VerSe CT benchmark (Sekuboyina et al., Med. Image Anal. 2021) underlies the
  VerSe meshes.
- **No redistribution.** The SpineToolkit repository has no licence file, so the files are
  **not redistributed** here. `scripts/fetch_data.py` downloads them and verifies each one.

## Files in this folder

| File | Content |
|---|---|
| `split_manifest_57.csv` | Vertebra-wise development split (20 calibration / 37 held-out) |
| `source_manifest.csv` | How each of the 1678 raw files is obtained (`download`, `copy`, `stl_from_vtk`) and its expected SHA-256 |
| `raw_file_sha256.csv` | SHA-256 of every raw input file used for the paper (layout `raw/VerSe_Ligaments_GT/...`) |
| `dissertation_archived_predictions/` | Predictions archived from the original analysis (provenance of the development tables) |

The subject codes used in every table are defined in `code/ligament/data_io.py` (`SUBJECT_CODE`).

## How the raw layout is rebuilt

`scripts/fetch_data.py` rebuilds the layout in three ways:

- **Downloaded unchanged (914 files).** The VTK meshes and annotations of S01–S03, and the STL
  meshes and annotations of S04.
- **Copied (712 files).** The `*_STL/*_ligaments*` folders. They hold copies of the official
  annotations, and all 712 files are byte-identical to the official ones. The loader reads
  these copies first.
- **Rebuilt as STL (52 files).** The STL surfaces of S01–S03 and S05, written from the official
  VTK meshes with one triangle per VTK face. The analyses read only their vertex arrays, and
  each array is verified against the files used for the paper. One of these STLs, `Verse581_T1`,
  differs from a plain conversion in the corner order of one degenerate triangle; the script
  reproduces that difference.

Only the development GCZ column (paper Table 6) and the original notebook read the STL files.
All other analyses use the VTK meshes.

## Notes

- **Frames.** Meshes and annotations of S01–S03 are in the LPS frame (VTK header `SPACE=LPS`;
  markup `coordinateSystem=LPS`). The S04 STL files declare `SPACE=RAS` while their annotations
  declare LPS. S04–S05 are predicted but not scored in the development analysis, and they are
  not part of the primary evaluation.
- **Short scans.** The S02 scan ends at L4.
