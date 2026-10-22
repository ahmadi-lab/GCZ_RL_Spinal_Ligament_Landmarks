"""Data loading for the ligament landmark study.

The loaders reproduce the file formats and search order of the original
analysis notebook (``fair_comparison.ipynb``):

* vertebral surface meshes: legacy binary VTK POLYDATA (``<pat>/<pat>_<lv>.vtk``)
  or STL (``<pat>_STL/<pat>_<lv>.stl``);
* expert landmarks: 3D Slicer markup files (``*.mrk.json``), one file per
  ligament sub-group (for example ``ALL_I`` and ``ALL_S``).

All coordinates are in millimetres in the LPS frame used by the source files.
"""
from __future__ import annotations

import json
import os
import re
from pathlib import Path

import numpy as np

PACKAGE_ROOT = Path(__file__).resolve().parents[2]
DATASET = PACKAGE_ROOT / "dataset"
# Raw meshes/annotations (scripts/fetch_data.py). GCZRL_RAW_ROOT may point elsewhere, e.g. to a
# short path on Windows, where the default location can exceed the 260-character path limit.
RAW_BASE = Path(os.environ.get("GCZRL_RAW_ROOT", DATASET / "raw"))
RAW_ROOT = RAW_BASE / "VerSe_Ligaments_GT"
MANIFEST_57 = DATASET / "split_manifest_57.csv"

LIGAMENTS = ["ALL", "PLL", "LF", "ISL", "SSL", "ITL", "CL"]
LIGAMENT_SUFFIXES = {
    "ALL": ["ALL_I", "ALL_S"],
    "PLL": ["PLL_I", "PLL_S"],
    "LF": ["LF_I", "LF_S"],
    "ISL": ["ISL_I", "ISL_S"],
    "SSL": ["SSL"],
    "ITL": ["ITL_L", "ITL_R"],
    "CL": ["CL_I_L", "CL_I_R", "CL_S_L", "CL_S_R"],
}
LEVEL_ORDER = ["C1", "C2", "C3", "C4", "C5", "C6", "C7",
               "T1", "T2", "T3", "T4", "T5", "T6", "T7", "T8", "T9", "T10", "T11", "T12", "T13",
               "L1", "L2", "L3", "L4", "L5"]

# Anonymized subject codes used in every reported table.  The mapping to the
# public VerSe / source identifiers is needed only to read the raw files.
SUBJECT_CODE = {"Verse005": "S01", "Verse561": "S02", "Verse581": "S03",
                "Bc11701": "S04", "Bc11702": "S05"}
SOURCE_ID = {v: k for k, v in SUBJECT_CODE.items()}


def read_vtk_polydata(path: Path) -> tuple[np.ndarray, np.ndarray]:
    """Read a legacy *binary* VTK POLYDATA file into (vertices, triangles)."""
    with open(path, "rb") as fh:
        while True:
            line = fh.readline()
            if not line:
                raise ValueError(f"No POINTS section in {path}")
            s = line.decode("ascii", errors="ignore").strip()
            if s.upper().startswith("POINTS"):
                n_pts = int(s.split()[1])
                break
        verts = np.frombuffer(fh.read(n_pts * 3 * 4), dtype=">f4").reshape(n_pts, 3).astype(float)
        n_vals = None
        while True:
            line = fh.readline()
            if not line:
                break
            s = line.decode("ascii", errors="ignore").strip()
            if s.upper().startswith("POLYGONS"):
                n_vals = int(s.split()[2])
                break
        raw = np.frombuffer(fh.read(n_vals * 4), dtype=">i4") if n_vals else np.zeros(0, ">i4")
    faces = []
    i = 0
    while i < len(raw):
        cnt = int(raw[i])
        if cnt == 3 and i + 3 < len(raw):
            faces.append(raw[i + 1:i + 4])
        i += cnt + 1
    return verts, np.array(faces, dtype=np.int64).reshape(-1, 3)


def read_stl(path: Path) -> tuple[np.ndarray, np.ndarray]:
    import trimesh
    mesh = trimesh.load(path, process=False)
    return np.asarray(mesh.vertices, dtype=float), np.asarray(mesh.faces, dtype=np.int64)


def mesh_path(patient: str, level: str, root: Path = RAW_ROOT) -> Path | None:
    """Mesh path using the original notebook search order (STL first, then VTK)."""
    stl = root / f"{patient}_STL" / f"{patient}_{level}.stl"
    if stl.exists():
        return stl
    vtk = root / patient / f"{patient}_{level}.vtk"
    if vtk.exists():
        return vtk
    return None


def load_mesh(patient: str, level: str, root: Path = RAW_ROOT):
    path = mesh_path(patient, level, root)
    if path is None:
        return None, None
    if path.suffix.lower() == ".stl":
        return read_stl(path)
    return read_vtk_polydata(path)


def landmark_dir(patient: str, level: str, root: Path = RAW_ROOT) -> Path | None:
    for sub in (f"{patient}_STL", patient):
        for suffix in ("_ligaments_GT", "_ligaments"):
            d = root / sub / f"{patient}_{level}{suffix}"
            if d.exists():
                return d
    return None


def read_markup(path: Path) -> np.ndarray | None:
    with open(path, encoding="utf-8") as fh:
        data = json.load(fh)
    pts = [cp["position"] for cp in data["markups"][0]["controlPoints"]
           if cp.get("positionStatus", "defined") == "defined"]
    return np.asarray(pts, dtype=float) if pts else None


def landmark_files(patient: str, level: str, root: Path = RAW_ROOT) -> dict[str, list[Path]]:
    d = landmark_dir(patient, level, root)
    out: dict[str, list[Path]] = {lig: [] for lig in LIGAMENTS}
    if d is None:
        return out
    names = sorted(p.name for p in d.iterdir())
    for lig, sfxs in LIGAMENT_SUFFIXES.items():
        for suf in sfxs:
            pat = re.compile(rf".*_{re.escape(suf)}(_\d+)?\.mrk\.json$")
            out[lig].extend(d / n for n in names if pat.match(n))
    return out


def load_landmarks(patient: str, level: str, root: Path = RAW_ROOT) -> dict | None:
    """Return {ligament: (N,3) array or None}, or None when no annotation folder exists."""
    if landmark_dir(patient, level, root) is None:
        return None
    result = {}
    for lig, files in landmark_files(patient, level, root).items():
        parts = [p for p in (read_markup(f) for f in files) if p is not None]
        result[lig] = np.vstack(parts) if parts else None
    return result


def read_manifest(path: Path = MANIFEST_57):
    import pandas as pd
    df = pd.read_csv(path)
    calib = [(r["patient"], r["level"]) for _, r in df.iterrows() if r["split"] == "calibration"]
    held = [(r["patient"], r["level"]) for _, r in df.iterrows() if r["split"] == "heldout"]
    return calib, held


def case_code(patient: str, level: str) -> str:
    return f"{SUBJECT_CODE.get(patient, patient)}_{level}"
