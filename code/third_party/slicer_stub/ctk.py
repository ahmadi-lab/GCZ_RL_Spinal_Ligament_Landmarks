"""Stub for Common Toolkit (Slicer-embedded)."""
class ctkCollapsibleButton:
    def __init__(self): self.text = ""
class ctkSliderWidget:
    value = 0.0
    def __init__(self): pass

# ============================================================
# STEP 11 NOTEBOOK CELL (copy/paste into GCZ_journal notebook)
# ============================================================
# Cell code to add after Step 10 in the main notebook:
#
# # ===========================================================
# # STEP 11 — Run Official SpineToolkit SLD on 37 Held-out STL
# # ===========================================================
# import importlib, sys, os
# import pandas as pd
#
# _WS = r"<workspace>"
# if _WS not in sys.path: sys.path.insert(0, _WS)
# import check_vtk as sld_runner
# importlib.reload(sld_runner)
#
# # -- Build STL file-path lookup ----------------------------
# # Adjust _STL_BASE to match where your per-vertebra STL files are stored.
# # They should be the same STL files used by the GCZ pipeline (Step 8).
# _STL_BASE = OUTPUT_BASE   # from Step 1 setup cell
#
# def _stl_fn(case_id, level):
#     """Return path to STL for case_id / level from the GCZ pipeline cache."""
#     # Try the mesh cache from Step 7.5 (variable _mesh75)
#     try:
#         return _mesh75[case_id][level]   # vtkPolyData -> write to temp STL
#     except Exception:
#         pass
#     # Fallback: look for STL file on disk
#     for ext in ['.stl', '.STL']:
#         p = os.path.join(_STL_BASE, case_id, 'segmentations', level + ext)
#         if os.path.exists(p): return p
#         p = os.path.join(_STL_BASE, case_id, level + ext)
#         if os.path.exists(p): return p
#     raise FileNotFoundError(f"STL not found for {case_id}/{level}")
#
# # -- Coordinate frame for VerSe/LPS data -------------------
# # VerSe STL files are in LPS (Left-Posterior-Superior).
# # SpineLib defaults to RAS. For LPS: Anterior = -y, so override:
# _SPINE_R = [1, 0, 0]   # Left  (same in LPS and RAS)
# _SPINE_A = [0, -1, 0]  # Anterior = -y in LPS
# _SPINE_S = [0, 0, 1]   # Superior (same in LPS and RAS)
#
# # -- Run SLD ------------------------------------------------
# _TEST_CASES = _TEST75   # list of (case_id, level_str) tuples from Step 7.5
#
# df_sld_official = sld_runner.run_all(
#     _TEST_CASES, _stl_fn,
#     spine_r=_SPINE_R, spine_a=_SPINE_A, spine_s=_SPINE_S,
#     max_angle=45.0
# )
#
# print(f"\nSLD-official predictions: {len(df_sld_official)} rows")
# print(df_sld_official.groupby('ligament').size())
#
# # Save
# _out = os.path.join(OUTPUT_BASE, 'sld_official_predictions_test.csv')
# df_sld_official.to_csv(_out, index=False)
# print(f"Saved: {os.path.basename(_out)}")
# ===========================================================

