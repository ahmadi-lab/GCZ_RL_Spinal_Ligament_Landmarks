"""Replacement for the 3D Slicer *ExtractCenterline* (SlicerVMTK) module.

The official SpineToolkit shape decomposition calls VMTK to obtain the medial
centerline of each vertebral process between two seed points.  VMTK is not
available for this Python environment (no pip wheel, no Slicer installation),
so this module provides a documented approximation:

1. the input surface is decimated to about ``targetNumberOfPoints`` vertices
   (as the Slicer module's preprocessing does);
2. a geodesic (Dijkstra) path is computed on the surface between the vertices
   closest to the two seeds;
3. every path point is moved to the centroid of the surface vertices lying in
   a thin slab perpendicular to the local path direction (cross-section
   centroid).  For tubular structures such as transverse, spinous and
   articular processes this centroid lies on the medial axis, which is what
   the VMTK centerline approximates;
4. the medial polyline is lightly smoothed and resampled.

This is not VMTK.  Results produced with it are labelled
"official SpineLib code, VMTK centerline replaced by a cross-section medial
approximation" throughout the package.  Set ``MEDIAL=False`` to reproduce the
earlier (June 2026) stub that used the raw surface geodesic path.
"""
import numpy as np
import vtk
from vtk.util.numpy_support import numpy_to_vtk, vtk_to_numpy

MEDIAL = True
SLAB_HALF_WIDTH_MM = 1.0
CROSS_SECTION_RADIUS_MM = 12.0
N_OUT = 40


def _resample(pts, n):
    pts = np.asarray(pts, float)
    if len(pts) < 2:
        return pts
    seg = np.linalg.norm(np.diff(pts, axis=0), axis=1)
    arc = np.concatenate([[0.0], np.cumsum(seg)])
    if arc[-1] < 1e-9:
        return np.repeat(pts[:1], n, axis=0)
    u = np.linspace(0, arc[-1], n)
    return np.column_stack([np.interp(u, arc, pts[:, i]) for i in range(3)])


def _medial(path, verts):
    path = _resample(path, N_OUT)
    tang = np.gradient(path, axis=0)
    norms = np.linalg.norm(tang, axis=1, keepdims=True)
    tang = tang / np.maximum(norms, 1e-9)
    out = path.copy()
    for i, (p, t) in enumerate(zip(path, tang)):
        rel = verts - p
        along = rel @ t
        radial = np.linalg.norm(rel - np.outer(along, t), axis=1)
        sel = (np.abs(along) <= SLAB_HALF_WIDTH_MM) & (radial <= CROSS_SECTION_RADIUS_MM)
        if sel.sum() >= 6:
            out[i] = verts[sel].mean(0)
    # keep the seeds as end points (the Slicer module seeds the centerline there)
    out[0], out[-1] = path[0], path[-1]
    sm = out.copy()
    for i in range(1, len(out) - 1):
        sm[i] = 0.25 * out[i - 1] + 0.5 * out[i] + 0.25 * out[i + 1]
    return sm


class ExtractCenterlineLogic:

    def preprocess(self, polydata, targetNumberOfPoints=5000.0, decimationAggressiveness=4,
                   subdivideInputSurface=False):
        n = polydata.GetNumberOfPoints()
        pd = polydata
        if n > targetNumberOfPoints * 1.2:
            dec = vtk.vtkQuadricDecimation()
            dec.SetInputData(polydata)
            dec.SetTargetReduction(float(min(0.95, 1.0 - targetNumberOfPoints / n)))
            dec.Update()
            pd = dec.GetOutput()
        clean = vtk.vtkCleanPolyData()
        clean.SetInputData(pd)
        clean.Update()
        tri = vtk.vtkTriangleFilter()
        tri.SetInputData(clean.GetOutput())
        tri.Update()
        return tri.GetOutput()

    def extractCenterline(self, preprocessedPolyData, pointMarkup):
        pts = getattr(pointMarkup, "control_points", None)
        if pts is None or len(pts) < 2:
            raise ValueError("centerline needs two seed points")
        start, end = np.asarray(pts[0], float), np.asarray(pts[1], float)
        loc = vtk.vtkPointLocator()
        loc.SetDataSet(preprocessedPolyData)
        loc.BuildLocator()
        a, b = loc.FindClosestPoint(start), loc.FindClosestPoint(end)
        path = None
        if a != b:
            dj = vtk.vtkDijkstraGraphGeodesicPath()
            dj.SetInputData(preprocessedPolyData)
            dj.SetStartVertex(a)
            dj.SetEndVertex(b)
            dj.Update()
            if dj.GetOutput().GetNumberOfPoints() >= 2:
                path = vtk_to_numpy(dj.GetOutput().GetPoints().GetData())[::-1]
        if path is None:
            path = np.linspace(start, end, 24)
        path = np.vstack([start, path, end])
        verts = vtk_to_numpy(preprocessedPolyData.GetPoints().GetData()).astype(float)
        line = _medial(path, verts) if MEDIAL else _resample(path, N_OUT)
        p = vtk.vtkPoints()
        p.SetData(numpy_to_vtk(np.ascontiguousarray(line), deep=True))
        out = vtk.vtkPolyData()
        out.SetPoints(p)
        return out, vtk.vtkPolyData()

    def createCurveTreeFromCenterline(self, centerlinePolyData, centerlineCurveNode, propertiesTableNode):
        if centerlinePolyData is None or centerlinePolyData.GetPoints() is None:
            return
        if centerlineCurveNode is not None:
            centerlineCurveNode.control_points = vtk_to_numpy(centerlinePolyData.GetPoints().GetData()).copy()
