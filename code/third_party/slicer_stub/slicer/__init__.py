"""
Minimal 3D Slicer stub for running SpineToolkit standalone (without Slicer installed).

Usage: add the slicer_stub/ directory to sys.path BEFORE importing SpineLib.
    import sys
    sys.path.insert(0, "/path/to/SpineToolkit/slicer_stub")
    import slicer, qt, ctk, ExtractCenterline
    import SpineLib
"""
import numpy as np
from scipy.interpolate import CubicSpline


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------
class _DisplayStub:
    def SetTextScale(self, *a): pass
    def SetSelectedColor(self, *a): pass
    def SetColor(self, *a): pass
    def SetOpacity(self, *a): pass
    def SetGlyphScale(self, *a): pass
    def SetVisibility(self, *a): pass
    def SetHandlesInteractive(self, *a): pass


def _resample_polyline(pts, n):
    pts = np.asarray(pts, dtype=float)
    if len(pts) < 2:
        return pts
    diffs = np.linalg.norm(np.diff(pts, axis=0), axis=1)
    arc = np.concatenate([[0.0], np.cumsum(diffs)])
    total = arc[-1]
    if total < 1e-9:
        return np.tile(pts[0], (n, 1))
    u_new = np.linspace(0, total, n)
    return np.column_stack([np.interp(u_new, arc, pts[:, i]) for i in range(3)])


def _resample_curve_polynomial(pts, n):
    pts = np.asarray(pts, dtype=float)
    if len(pts) < 2:
        return pts
    if len(pts) == 2:
        return np.linspace(pts[0], pts[1], n)
    diffs = np.linalg.norm(np.diff(pts, axis=0), axis=1)
    arc = np.concatenate([[0.0], np.cumsum(diffs)])
    total = arc[-1]
    if total < 1e-9:
        return np.tile(pts[0], (n, 1))
    _, keep = np.unique(arc, return_index=True)
    arc_u = arc[keep] / total
    pts_u = pts[keep]
    try:
        cs = CubicSpline(arc_u, pts_u, bc_type="not-a-knot")
        return cs(np.linspace(0, 1, n))
    except Exception:
        return _resample_polyline(pts, n)


class _PointsAccessor:
    def __init__(self, pts):
        self._pts = np.asarray(pts, dtype=float)
    def GetData(self):
        # vtk_to_numpy expects a VTK array; wrap numpy array
        from vtk.util.numpy_support import numpy_to_vtk
        if self._pts.size == 0:
            import vtk as _vtk
            arr = _vtk.vtkFloatArray()
            arr.SetNumberOfComponents(3)
            return arr
        return numpy_to_vtk(self._pts.astype(np.float64), deep=True)


class _NodeStub:
    """Generic no-op Slicer scene node."""
    def __init__(self, class_name="", name=""):
        self._class = class_name
        self.name = name
        self._polydata = None
        self.control_points = np.zeros((0, 3))

    def _add_pt(self, pt):
        pt = np.asarray(pt, dtype=float).reshape(1, 3)
        self.control_points = (np.vstack([self.control_points, pt])
                               if self.control_points.shape[0] else pt)

    def AddControlPoint(self, pt): self._add_pt(pt)
    def AddFiducialFromArray(self, pt, *a): self._add_pt(pt)
    def AddControlPointWorld(self, pt, *a): self._add_pt(pt)
    def GetDisplayNode(self): return _DisplayStub()
    def SetAndObservePolyData(self, pd): self._polydata = pd
    def SetAndObserveTransformNodeID(self, *a): pass
    def HardenTransform(self): pass
    def SetMatrixTransformToParent(self, *a): pass
    def SetOrigin(self, *a): pass
    def SetNormal(self, *a): pass
    def SetSize(self, *a): pass
    def SetCenter(self, *a): pass
    def GetID(self): return ""
    def SetCurveTypeToPolynomial(self): pass
    def SetCurveTypeToCardinalSpline(self): pass
    def GetCurveLengthWorld(self):
        # Slicer returns the arc length of the curve; the earlier stub returned a
        # constant 1.0, which made resampleCurve() emit hundreds of points.
        pts = np.asarray(self.control_points, dtype=float)
        if pts.shape[0] < 2:
            return 0.0
        return float(np.linalg.norm(np.diff(pts, axis=0), axis=1).sum())
    def GetCurvePointIndexFromControlPointIndex(self, i): return i
    def GetCurvePointToWorldTransformAtPointIndex(self, i, mat): pass

    def ResampleCurveWorld(self, sample_dist):
        pts = self.control_points
        if pts.shape[0] < 2 or sample_dist < 1e-9:
            return
        diffs = np.linalg.norm(np.diff(pts, axis=0), axis=1)
        total = float(diffs.sum())
        if total < 1e-9:
            return
        n_pts = max(2, int(round(total / sample_dist)) + 1)
        self.control_points = _resample_polyline(pts, n_pts)

    def GetCurvePointsWorld(self):
        return _PointsAccessor(self.control_points)

    def GetPolyData(self):
        return self._polydata


# ---------------------------------------------------------------------------
# slicer.mrmlScene
# ---------------------------------------------------------------------------
class _MRMLScene:
    def AddNewNodeByClass(self, class_name, name=""):
        return _NodeStub(class_name, name)
    def RemoveNode(self, node): pass
    def GetNodesByClass(self, class_name): return []

mrmlScene = _MRMLScene()


# ---------------------------------------------------------------------------
# slicer.app
# ---------------------------------------------------------------------------
class _App:
    class _LayoutMgr:
        def setLayout(self, *a): pass
    def processEvents(self): pass
    def layoutManager(self): return self._LayoutMgr()

app = _App()


# ---------------------------------------------------------------------------
# slicer.util
# ---------------------------------------------------------------------------
class util:
    @staticmethod
    def arrayFromMarkupsControlPoints(node):
        if node is None:
            return np.zeros((0, 3))
        pts = getattr(node, "control_points", None)
        if pts is None or (hasattr(pts, "__len__") and len(pts) == 0):
            return np.zeros((0, 3))
        return np.asarray(pts, dtype=float)

    @staticmethod
    def updateMarkupsControlPointsFromArray(node, pts):
        if node is not None:
            node.control_points = np.asarray(pts, dtype=float)

    @staticmethod
    def arrayFromModelPoints(node):
        if node is not None and getattr(node, "_polydata", None) is not None:
            from vtk.util.numpy_support import vtk_to_numpy
            return vtk_to_numpy(node._polydata.GetPoints().GetData())
        return np.zeros((0, 3))

    @staticmethod
    def getNodesByClass(class_name): return []
    @staticmethod
    def getNode(name): return None
    @staticmethod
    def loadModel(path): return _NodeStub("vtkMRMLModelNode")
    @staticmethod
    def loadColorTable(path): pass
    @staticmethod
    def plot(*a, **kw): pass


# ---------------------------------------------------------------------------
# Standalone createResampledCurve — monkey-patched into SpineLib.SlicerTools
# ---------------------------------------------------------------------------
def _standalone_createResampledCurve(control_points, numberOfPoints,
                                     name="MarkupsCurve", color=None):
    pts = np.asarray(control_points, dtype=float)
    resampled = _resample_curve_polynomial(pts, numberOfPoints)
    node = _NodeStub("vtkMRMLMarkupsCurveNode", name)
    node.control_points = resampled
    return node


# ---------------------------------------------------------------------------
# ScriptedLoadableModule stubs (only needed if LigamentDetection.py is imported)
# ---------------------------------------------------------------------------
class ScriptedLoadableModule:
    def __init__(self, parent): pass
    def getDefaultModuleDocumentationLink(self): return ""

class ScriptedLoadableModuleWidget:
    layout = None
    def setup(self): pass
    def cleanup(self): pass

class ScriptedLoadableModuleLogic:
    pass

class ScriptedLoadableModuleTest:
    pass


# pip_install and extensionsManagerModel — no-ops in standalone mode
import subprocess, sys as _sys

class _ExtManagerModel:
    interactive = False
    def installExtensionFromServer(self, name, restart=True): return True

# Add methods to existing classes via monkey-patch after class definition
def _pip_install_stub(package_or_req):
    """No-op: packages already installed via requirements.txt."""
    pass

# Patch util
util.pip_install = staticmethod(_pip_install_stub)

# Patch app
app.extensionsManagerModel = lambda: _ExtManagerModel()


# ---------------------------------------------------------------------------
# Catch-all for slicer.vtkMRML* and other unrecognized attributes
# ---------------------------------------------------------------------------
def __getattr__(name):
    """Return a generic stub class for any unrecognized slicer.* attribute."""
    return type(name, (), {"__init__": lambda self, *a, **kw: None})
