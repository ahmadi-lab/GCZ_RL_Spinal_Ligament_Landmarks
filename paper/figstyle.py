"""Shared style for the figures of the paper: typography, palette, markers, and the 3D renderer.

Methods are distinguished by colour and marker shape (GCZ square, GCZ-RL circle, SLD triangle),
ligament classes by the Okabe-Ito palette. Quantitative figures are written as PDF and SVG with
text kept as text; anatomical renders are raster images placed in a vector page.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

LIGAMENTS = ["ALL", "PLL", "LF", "ISL", "SSL", "ITL", "CL"]
CLASS = {"ALL": "#0072B2", "PLL": "#56B4E9", "LF": "#009E73", "ISL": "#E69F00",
         "SSL": "#D55E00", "ITL": "#CC79A7", "CL": "#7B3294"}
METHOD = {"GCZ": "#3B6EA8", "GCZ-RL": "#C4513A", "SLD": "#4E9A6B"}
MARKER = {"GCZ": "s", "GCZ-RL": "o", "SLD": "^"}
INK, MUTED, GRID, NEUTRAL = "#1f1f1f", "#5f5f5f", "#e6e6e6", "#a6a6a6"
BONE, EXPERT = "#e8e2d6", "#2b2b2b"
# LPS frame: +x left, +y posterior, +z superior; (view direction, view-up)
VIEWS = {"superior": ((0, 0, 1), (0, -1, 0)), "left": ((1, 0, 0), (0, 0, 1)),
         "posterior": ((0, 1, 0), (0, 0, 1))}
VIEW_TITLE = {"superior": "Superior view", "left": "Left lateral view", "posterior": "Posterior view"}


def setup():
    plt.rcParams.update({
        "font.family": "Arial", "font.size": 7.5, "axes.titlesize": 8, "axes.labelsize": 7.5,
        "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7,
        "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
        "xtick.major.size": 2.5, "ytick.major.size": 2.5, "axes.edgecolor": INK,
        "axes.labelcolor": INK, "xtick.color": INK, "ytick.color": INK, "text.color": INK,
        "svg.fonttype": "none", "pdf.fonttype": 42, "savefig.dpi": 600, "figure.dpi": 150,
        "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False})


def panel(ax, letter, title="", x=-0.02):
    ax.set_title(f"({letter}) {title}".rstrip(), loc="left", fontsize=8, fontweight="bold", x=x, pad=4)


def ygrid(ax):
    ax.yaxis.grid(True, color=GRID, linewidth=0.5)
    ax.set_axisbelow(True)


def save(fig, out: Path, stem: str, vector=True):
    out.mkdir(parents=True, exist_ok=True)
    exts = ("pdf", "svg", "png") if vector else ("pdf", "png")
    for ext in exts:
        fig.savefig(out / f"{stem}.{ext}", dpi=600, bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)
    return [out / f"{stem}.{e}" for e in exts]


# ------------------------------------------------------------------------------ 3D renderer
def render(verts, faces, view, layers=(), zoom=1.25, focal=None, size=(1400, 1400),
           vertex_rgb=None, crop=True):
    """Off-screen render of a vertebral mesh with point, glyph, line, and mesh layers.

    layers: dicts with kind in {"spheres", "cubes", "lines", "mesh"}; see figures_main.py.
    """
    import pyvista as pv
    mesh = pv.PolyData(np.asarray(verts, float), np.hstack([np.full((len(faces), 1), 3), faces]).ravel())
    pl = pv.Plotter(off_screen=True, window_size=size)
    pl.set_background("white")
    if vertex_rgb is not None:
        mesh.point_data["rgb"] = vertex_rgb
        pl.add_mesh(mesh, scalars="rgb", rgb=True, smooth_shading=True, specular=0.05, ambient=0.25)
    else:
        pl.add_mesh(mesh, color=BONE, smooth_shading=True, specular=0.05, ambient=0.25)
    for lay in layers:
        kind = lay["kind"]
        pts = np.atleast_2d(np.asarray(lay.get("points", np.zeros((0, 3))), float))
        if kind == "spheres" and len(pts):
            pl.add_points(pts, color=lay["color"], point_size=lay.get("size", 14), render_points_as_spheres=True)
        elif kind == "cubes" and len(pts):
            s = lay.get("edge", 1.8)
            pl.add_mesh(pv.PolyData(pts).glyph(geom=pv.Cube(x_length=s, y_length=s, z_length=s), scale=False,
                                               orient=False), color=lay["color"], smooth_shading=False,
                        opacity=lay.get("opacity", 1.0))
        elif kind == "lines":
            for a, b in lay["pairs"]:
                if np.linalg.norm(np.asarray(b) - np.asarray(a)) > 0.05:
                    pl.add_mesh(pv.Line(a, b).tube(radius=lay.get("radius", 0.18)), color=lay["color"])
        elif kind == "mesh":
            m = pv.PolyData(np.asarray(lay["verts"], float),
                            np.hstack([np.full((len(lay["faces"]), 1), 3), lay["faces"]]).ravel())
            pl.add_mesh(m, color=lay["color"], opacity=lay.get("opacity", 0.4), smooth_shading=True)
    vec, up = VIEWS[view]
    pl.view_vector(vec, viewup=up)
    if focal is not None:
        pl.camera.focal_point = tuple(np.asarray(focal, float))
    pl.camera.zoom(zoom)
    img = pl.screenshot(return_img=True)
    pl.close()
    if crop:
        rows, cols = np.where(img.min(axis=2) < 250)
        if len(rows):
            img = img[max(rows.min() - 8, 0):rows.max() + 8, max(cols.min() - 8, 0):cols.max() + 8]
    return img


def region_colors(verts, regions: dict, order=("ALL", "CL", "LF", "PLL", "ISL", "SSL", "ITL"), mix=0.75):
    """Per-vertex RGB: bone colour, tinted by the class of each region (later classes on top)."""
    from matplotlib.colors import to_rgb
    base = np.tile(np.array(to_rgb(BONE)), (len(verts), 1))
    for lig in order:
        mask = regions.get(lig)
        if mask is not None:
            base[mask] = (1 - mix) * base[mask] + mix * np.array(to_rgb(CLASS[lig]))
    return (base * 255).astype(np.uint8)


def show(ax, img, title=None):
    ax.imshow(img)
    ax.set_anchor("N")        # images aligned at the top, so panel titles share one baseline
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    if title:
        ax.set_title(title, fontsize=7.5, pad=2)


def legend_handles(items):
    """items: list of (label, color, marker, size)."""
    from matplotlib.lines import Line2D
    return [Line2D([], [], linestyle="", marker=m, color=c, markersize=s, label=l, markeredgewidth=0)
            for l, c, m, s in items]
