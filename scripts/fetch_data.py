"""Download and rebuild the raw meshes and annotations used in the paper.

The meshes and expert ligament annotations are those published by the SLD authors in the
SpineToolkit repository (https://github.com/VisSim-UniKO/SpineToolkit, folder
LigamentLandmarkDetection/Datasets). That repository has no licence file, so this project does
not redistribute the files; this script fetches them from a pinned commit and rebuilds the exact
input layout used by the analyses:

  download      files fetched unchanged from the pinned commit
  copy          byte-identical copies of downloaded annotation files (the *_STL folders)
  stl_from_vtk  binary STL surfaces written from the official VTK meshes, one triangle per VTK
                face; the analyses read their vertex arrays, which are verified against
                the hashes of the files used for the paper

Every file is checked against dataset/source_manifest.csv (SHA-256 of the file, or of the
float64 vertex array for rebuilt STL files).  Exit code 1 if anything does not match.

Usage:
    python scripts/fetch_data.py                 # into dataset/raw
    python scripts/fetch_data.py --out D:/gczrl_raw   # elsewhere; then set GCZRL_RAW_ROOT
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import shutil
import struct
import sys
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "code" / "ligament"))
from data_io import read_stl, read_vtk_polydata  # noqa: E402

COMMIT = "62f39cffdc3957f2de35b5647c4e7bc9c79f9d28"
BASE_URL = f"https://raw.githubusercontent.com/VisSim-UniKO/SpineToolkit/{COMMIT}/LigamentLandmarkDetection/Datasets/"
MANIFEST = REPO / "dataset" / "source_manifest.csv"
# The STL used for the paper differs from a plain conversion in the corner order of one
# degenerate triangle; reproduce it exactly.
FACE_ORDER_PATCH = {"VerSe_Ligaments_GT/Verse581_STL/Verse581_T1.stl": {38747: (2, 1, 0)}}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def vertex_hash(path: Path) -> str:
    v, _ = read_stl(path)
    return hashlib.sha256(np.ascontiguousarray(v, dtype="<f8").tobytes()).hexdigest()


def download(src: str, dst: Path, retries: int = 3) -> None:
    dst.parent.mkdir(parents=True, exist_ok=True)
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(BASE_URL + src, timeout=120) as r, open(dst, "wb") as fh:
                shutil.copyfileobj(r, fh)
            return
        except OSError:
            if attempt == retries - 1:
                raise


def write_stl(vtk_path: Path, dst: Path, patch: dict) -> None:
    verts, faces = read_vtk_polydata(vtk_path)
    tri = verts[faces].astype(np.float32)
    for i, order in patch.items():
        tri[i] = tri[i][list(order)]
    n = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
    norm = np.linalg.norm(n, axis=1, keepdims=True)
    n = np.divide(n, norm, out=np.zeros_like(n), where=norm > 0).astype(np.float32)
    rec = np.zeros(len(tri), dtype=[("n", "<f4", 3), ("v", "<f4", (3, 3)), ("a", "<u2")])
    rec["n"], rec["v"] = n, tri
    dst.parent.mkdir(parents=True, exist_ok=True)
    with open(dst, "wb") as fh:
        fh.write(b"STL rebuilt by GCZ-RL scripts/fetch_data.py from the SpineToolkit VTK mesh".ljust(80, b" "))
        fh.write(struct.pack("<I", len(tri)))
        fh.write(rec.tobytes())


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, default=REPO / "dataset" / "raw")
    ap.add_argument("--workers", type=int, default=8)
    args = ap.parse_args()
    out = args.out
    rows = list(csv.DictReader(open(MANIFEST, encoding="utf-8")))
    sources = sorted({r["source"] for r in rows if r["kind"] in ("download", "copy", "stl_from_vtk")})
    cache = out / "_downloads"
    todo = [s for s in sources if not (cache / s).exists()]
    print(f"downloading {len(todo)} of {len(sources)} source files from SpineToolkit@{COMMIT[:7]} ...", flush=True)
    with ThreadPoolExecutor(args.workers) as ex:
        list(ex.map(lambda s: download(s, cache / s), todo))

    bad = []
    for r in rows:
        dst, src = out / r["target"], cache / r["source"]
        if r["kind"] in ("download", "copy"):
            dst.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(src, dst)
            ok = sha256(dst) == r["sha256"]
        else:
            write_stl(src, dst, FACE_ORDER_PATCH.get(r["target"], {}))
            ok = vertex_hash(dst) == r["sha256"]
        if not ok:
            bad.append(r["target"])
    shutil.rmtree(cache, ignore_errors=True)
    print(f"{len(rows) - len(bad)}/{len(rows)} files verified in {out}")
    if bad:
        print("MISMATCH:", *bad[:20], sep="\n  ")
        return 1
    if out.resolve() != (REPO / "dataset" / "raw").resolve():
        print(f"set GCZRL_RAW_ROOT={out} before running the analyses")
    return 0


if __name__ == "__main__":
    sys.exit(main())
