"""
vtk compat shim: inject vtkIdFilter alias into the vtk module
so that vtk_convenience.py (which imports it but never uses it) can load.
Run: import vtk_compat_shim  before  import vtk_convenience
"""
import vtk as _vtk

if not hasattr(_vtk, 'vtkIdFilter'):
    # vtkIdFilter was renamed/removed in VTK 9.x.
    # It is only imported (never called) in vtk_convenience.py, so any class works.
    class _vtkIdFilterCompat:
        """Compatibility stub for vtkIdFilter (not available in VTK 9.x)."""
        pass
    _vtk.vtkIdFilter = _vtkIdFilterCompat
    import vtkmodules.vtkFiltersCore as _fc
    _fc.vtkIdFilter = _vtkIdFilterCompat

# vtkParametricSpline and vtkParametricFunctionSource may not be exposed
# at the top-level vtk namespace in all VTK 9.x builds.
if not hasattr(_vtk, 'vtkParametricSpline'):
    try:
        from vtk import vtkParametricSpline as _vPS
        _vtk.vtkParametricSpline = _vPS
    except ImportError:
        from vtkmodules.vtkCommonComputationalGeometry import vtkParametricSpline as _vPS
        _vtk.vtkParametricSpline = _vPS

if not hasattr(_vtk, 'vtkParametricFunctionSource'):
    try:
        from vtk import vtkParametricFunctionSource as _vPFS
        _vtk.vtkParametricFunctionSource = _vPFS
    except ImportError:
        from vtkmodules.vtkFiltersSources import vtkParametricFunctionSource as _vPFS
        _vtk.vtkParametricFunctionSource = _vPFS
