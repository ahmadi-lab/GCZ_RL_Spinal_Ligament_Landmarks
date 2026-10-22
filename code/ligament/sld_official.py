"""Run the official SpineToolkit ligament landmark detection (SLD) outside 3D Slicer.

What is official
    * ``third_party/SpineLib`` is an unmodified copy of SpineLib from the public
      SpineToolkit repository (Blomenkamp et al., University of Koblenz).
    * Preprocessing follows ``LigamentDetection.py`` exactly: pyACVD remeshing
      (subdivide=3, clusters=5000), 10 smoothing iterations, hole filling
      (max 1000), cleaning; ``Spine(max_angle=60)``; shape decomposition with
      ``with_lamina=False``; ``LigamentLandmarks``; final snapping of every
      landmark to the closest point of the *input* model.

What is substituted (documented, not tuned)
    * Slicer scene/markup objects are replaced by a minimal stub
      (``third_party/slicer_stub``).
    * The VMTK centerline is replaced by a cross-section medial approximation
      (``third_party/slicer_stub/ExtractCenterline.py``), because VMTK is not
      installable in this environment.

Coordinate handling (the correction relative to the earlier standalone runs)
    The source meshes and annotations are stored in LPS (VTK header
    ``SPACE=LPS``; markup ``coordinateSystem=LPS``).  3D Slicer converts LPS
    model files to its internal RAS frame on loading (x -> -x, y -> -y), and
    SpineLib hard-codes the RAS axes (R=+x, A=+y, S=+z) when orienting each
    vertebra.  Feeding LPS coordinates directly therefore swaps anterior with
    posterior and left with right.  Here every mesh is converted to RAS before
    SLD runs, and every output landmark is converted back to LPS for
    evaluation.  No parameter is fitted to annotations.
"""
from __future__ import annotations

import contextlib
import io
import os
import sys
import time
import traceback
from pathlib import Path

import numpy as np

CODE = Path(__file__).resolve().parents[1]
STUB = CODE / "third_party" / "slicer_stub"
LIB_PARENT = CODE / "third_party"

LIB_IDS = ["L5", "L4", "L3", "L2", "L1", "T13", "T12", "T11", "T10", "T9", "T8", "T7", "T6",
           "T5", "T4", "T3", "T2", "T1", "C7", "C6", "C5", "C4", "C3", "C2", "C1"]
LPS_TO_RAS = np.array([-1.0, -1.0, 1.0])

KEY_MAP = {
    "ALL": ["inferiorALL", "superiorALL"],
    "PLL": ["inferiorPLL", "superiorPLL"],
    "LF": ["LF_sup", "LF_inf"],
    "ISL": ["ISL"],
    "SSL": ["SSL"],
    "ITL": ["ITL"],
    "CL": ["CL_left_sup", "CL_right_sup", "CL_left_inf", "CL_right_inf"],
}

_LOADED = None


def _load_spinelib(medial: bool = True):
    global _LOADED
    if _LOADED is not None:
        _LOADED[2].MEDIAL = medial      # the centerline variant must follow every call
        return _LOADED
    os.environ.setdefault("MPLBACKEND", "Agg")
    os.environ.setdefault("matplotlibback", "Agg")
    import matplotlib
    matplotlib.use("Agg")
    for p in (str(STUB), str(LIB_PARENT)):
        if p not in sys.path:
            sys.path.insert(0, p)
    import vtk_compat_shim  # noqa: F401  (vtkIdFilter alias for VTK 9)
    import ExtractCenterline
    ExtractCenterline.MEDIAL = medial
    with contextlib.redirect_stdout(io.StringIO()):
        import SpineLib
    import vtk_convenience as conv
    _LOADED = (SpineLib, conv, ExtractCenterline)
    return _LOADED


def _polydata(verts, faces):
    import vtk
    from vtk.util.numpy_support import numpy_to_vtk, numpy_to_vtkIdTypeArray
    pts = vtk.vtkPoints()
    pts.SetData(numpy_to_vtk(np.ascontiguousarray(verts, dtype=np.float64), deep=True))
    cells = np.hstack([np.full((len(faces), 1), 3, dtype=np.int64), faces.astype(np.int64)]).ravel()
    ca = vtk.vtkCellArray()
    ca.SetCells(len(faces), numpy_to_vtkIdTypeArray(cells, deep=True))
    pd = vtk.vtkPolyData()
    pd.SetPoints(pts)
    pd.SetPolys(ca)
    return pd


def _as_points(x):
    if x is None:
        return None
    try:
        arr = np.asarray(x, dtype=float)
        if arr.ndim == 1 and arr.size == 3:
            arr = arr.reshape(1, 3)
        if arr.ndim == 2 and arr.shape[1] == 3 and len(arr):
            return arr
    except (ValueError, TypeError):
        parts = []
        for item in x:
            try:
                parts.append(np.asarray(item, dtype=float).reshape(-1, 3))
            except Exception:
                pass
        if parts:
            return np.vstack(parts)
    return None


def run_subject(levels_meshes: dict, medial: bool = True, log=print) -> dict:
    """Run SLD on the vertebrae of one subject.

    ``levels_meshes`` maps level -> (vertices LPS, faces).  Returns
    {level: {ligament: (N,3) LPS or None}} plus diagnostics under ``"_diag"``.
    """
    SpineLib, conv, _ = _load_spinelib(medial)
    import matplotlib.pyplot as plt
    t0 = time.time()
    levels = sorted(levels_meshes, key=LIB_IDS.index)          # inferior -> superior
    raw_ras = {lv: _polydata(levels_meshes[lv][0] * LPS_TO_RAS, levels_meshes[lv][1]) for lv in levels}
    geoms = []
    diag = {"errors": [], "levels": levels}
    for lv in levels:
        g = conv.polydata_remesh(raw_ras[lv], subdivide=3, clusters=5000)
        g = conv.polydata_smooth(g, iterations=10)
        g = conv.polydata_fillHoles(g, maximumHoleSize=1000.0)
        g = conv.polydata_clean(g)
        geoms.append(g)
    pbm = SpineLib.ProgressBarManager(layout=None)   # official class; no dialog is created
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        spine = SpineLib.Spine(geometries=geoms, indices=[LIB_IDS.index(lv) for lv in levels], max_angle=60.0)
        for vt in spine.vertebrae:
            try:
                plt.close("all")
                vt.get_shape_decomposition(progressBarManager=pbm, with_lamina=False)
            except Exception as e:  # official module also continues on failure
                diag["errors"].append(f"shape_decomposition {vt.name}: {type(e).__name__}: {e}")
        try:
            plt.close("all")
            SpineLib.LigamentLandmarks(spine, progressBarManager=pbm)
        except Exception as e:
            diag["errors"].append(f"LigamentLandmarks: {type(e).__name__}: {e}")
    diag["stdout_tail"] = buf.getvalue()[-4000:]
    out = {}
    for vt in spine.vertebrae:
        lv = vt.name
        lm = vt.ligament_landmarks or {}
        diag.setdefault("keys", {})[lv] = sorted(lm)
        res = {}
        for lig, keys in KEY_MAP.items():
            parts = []
            for k in keys:
                pts = _as_points(lm.get(k))
                if pts is None:
                    continue
                try:  # as in LigamentDetection.py: snap to the closest point of the input model
                    pts = np.asarray(conv.find_closest_points(raw_ras[lv], pts), float).reshape(-1, 3)
                except Exception as e:
                    diag["errors"].append(f"snap {lv} {k}: {e}")
                parts.append(pts)
            res[lig] = np.vstack(parts) * LPS_TO_RAS if parts else None
        out[lv] = res
    diag["seconds"] = round(time.time() - t0, 1)
    diag["centerline_variant"] = "medial" if _LOADED[2].MEDIAL else "surface_geodesic"
    out["_diag"] = diag
    return out
