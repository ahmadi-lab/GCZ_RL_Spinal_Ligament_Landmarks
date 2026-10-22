"""Cached loading of meshes and annotations for the 57 vertebrae."""
from __future__ import annotations

from functools import lru_cache

import numpy as np

from data_io import RAW_ROOT, load_landmarks, read_stl, read_vtk_polydata


@lru_cache(maxsize=None)
def mesh(patient: str, level: str, source: str = "stl") -> np.ndarray | None:
    """Vertices of one vertebra.

    ``source='stl'`` reproduces notebook cell 3 (STL first, VTK fallback), used by
    the original held-out GCZ run.  ``source='vtk'`` reproduces the later GCZ-RL
    cells, which read the SpineToolkit VTK files first.  Both files describe the
    same surface; the STL stores each vertex once per incident face.
    """
    stl = RAW_ROOT / f"{patient}_STL" / f"{patient}_{level}.stl"
    vtk = RAW_ROOT / patient / f"{patient}_{level}.vtk"
    order = [stl, vtk] if source == "stl" else [vtk, stl]
    for path in order:
        if path.exists():
            v = read_stl(path)[0] if path.suffix == ".stl" else read_vtk_polydata(path)[0]
            v.setflags(write=False)
            return v
    return None


@lru_cache(maxsize=None)
def landmarks(patient: str, level: str):
    return load_landmarks(patient, level)


def train_dict(keys, source: str) -> dict:
    return {k: {"verts": mesh(*k, source), "gt": landmarks(*k)} for k in keys if mesh(*k, source) is not None}
