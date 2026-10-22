"""Stub for slicer.util — all no-ops for standalone SpineLib use."""
import numpy as np


def _get_node_pts(node):
    if node is None:
        return np.zeros((0, 3))
    pts = getattr(node, "control_points", None)
    if pts is None or (hasattr(pts, "__len__") and len(pts) == 0):
        return np.zeros((0, 3))
    return np.asarray(pts, dtype=float)


def arrayFromMarkupsControlPoints(node):
    return _get_node_pts(node)

def updateMarkupsControlPointsFromArray(node, pts):
    if node is not None:
        node.control_points = np.asarray(pts, dtype=float)

def arrayFromModelPoints(node):
    if node is not None and getattr(node, "_polydata", None) is not None:
        from vtk.util.numpy_support import vtk_to_numpy
        return vtk_to_numpy(node._polydata.GetPoints().GetData())
    return np.zeros((0, 3))

def pip_install(package_or_req):
    """No-op: packages already installed."""
    pass

def getNodesByClass(class_name): return []
def getNode(name): return None
def loadModel(path):
    from slicer import _NodeStub
    return _NodeStub("vtkMRMLModelNode")
def loadColorTable(path): pass
def plot(*a, **kw): pass
