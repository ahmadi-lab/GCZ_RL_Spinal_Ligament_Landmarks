"""GCZ-RL: Geometry-Calibrated Zoning with bounded Residual Learning (Extra Trees).

For every GCZ landmark point p_GCZ of ligament class L on vertebra v:

    phi        = point_features(v, p_GCZ, L)              36 geometric features (no annotations)
    dp         = F_L(phi)                                  Extra Trees regressor (3-D output)
    dp         = clip(dp, tau_L)                           |dp| <= tau_L (correction radius)
    p_corr     = p_GCZ + alpha_L * dp_projected            blend weight alpha_L
    p_final    = nearest vertex of v inside the calibrated zone of L

Training data (from training vertebrae only)
    For each training vertebra, GCZ is run without that vertebra's own annotations; each
    GCZ point is paired with an expert landmark by Hungarian assignment and the difference
    (expert - GCZ) is the residual target.  Bilateral ligaments (CL, ITL, LF, ISL) are
    augmented by mirroring point and residual across the vertebral midline.
    tau_L = mean + 2 SD of the training residual magnitudes of class L.

Features (36): AABB-normalised position of the GCZ point, distance to the midline,
distance to the zone centroid and boundary, position inside the zone, zone centroid,
spread, size and principal-axis shares, vertebral level, ligament one-hot (7),
local surface normal (3), curvature, signed side, and the distance and direction to the
training-set mean landmark position of class L ("benchmark prior", 4).
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.optimize import linear_sum_assignment
from scipy.spatial import cKDTree
from sklearn.ensemble import (ExtraTreesRegressor, GradientBoostingRegressor,
                              HistGradientBoostingRegressor, RandomForestRegressor)
from sklearn.linear_model import Ridge
from sklearn.multioutput import MultiOutputRegressor
from sklearn.neighbors import KNeighborsRegressor
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

import gcz
from data_io import LEVEL_ORDER, LIGAMENTS

BILATERAL = {"CL", "ITL", "LF", "ISL"}
LVL_IDX = {l: i for i, l in enumerate(LEVEL_ORDER)}
N_JOBS = -1

# --------------------------------------------------------------------------- regressors
EXTRA_TREES = {
    "ET-1": dict(n_estimators=220, max_depth=18, min_samples_leaf=1, max_features=0.85,
                 bootstrap=False, random_state=8),
    "ET-2": dict(n_estimators=520, max_depth=None, min_samples_leaf=1, max_features=1.0,
                 bootstrap=False, random_state=91),
    "ET-3": dict(n_estimators=520, max_depth=13, min_samples_leaf=2, max_features=0.75,
                 bootstrap=False, random_state=92),
    "ET-4": dict(n_estimators=520, max_depth=18, min_samples_leaf=1, max_features=0.85,
                 bootstrap=True, random_state=93),
}


def make_regressor(name):
    """All candidate residual regressors used in the dissertation (Table 4.3)."""
    if name in EXTRA_TREES:
        return ExtraTreesRegressor(n_jobs=N_JOBS, **EXTRA_TREES[name])
    if name == "Ridge":
        return make_pipeline(StandardScaler(), Ridge(alpha=25.0))
    if name == "KNN":
        return "knn"  # neighbour count depends on the training size; see fit()
    if name == "RandomForest":
        return RandomForestRegressor(n_estimators=420, max_depth=None, min_samples_leaf=1, max_features=0.85,
                                     bootstrap=True, random_state=94, n_jobs=N_JOBS)
    if name == "GradientBoosting":
        return MultiOutputRegressor(GradientBoostingRegressor(
            n_estimators=260, learning_rate=0.035, max_depth=3, subsample=0.85, min_samples_leaf=2,
            random_state=95), n_jobs=N_JOBS)
    if name == "HistGradientBoosting":
        # one output at a time and one OpenMP thread (see fit): the multithreaded histogram sums of
        # HistGradientBoosting change in the last digits from run to run
        return make_pipeline(StandardScaler(), MultiOutputRegressor(HistGradientBoostingRegressor(
            max_iter=260, learning_rate=0.035, max_leaf_nodes=31, l2_regularization=0.05, min_samples_leaf=8,
            random_state=96), n_jobs=1))
    raise ValueError(name)


def fit(name, X, Y):
    model = make_regressor(name)
    if model == "knn":
        k = min(31, max(5, int(np.sqrt(len(X)))))
        model = make_pipeline(StandardScaler(), KNeighborsRegressor(n_neighbors=k, weights="distance"))
    if name == "HistGradientBoosting":
        from threadpoolctl import threadpool_limits
        with threadpool_limits(limits=1, user_api="openmp"):
            return model.fit(X, Y)
    return model.fit(X, Y)


# ----------------------------------------------------------------- final configuration
# Final GCZ-RL configuration reported in the dissertation (Chapter 4, Table 4.4):
# one Extra Trees regressor per ligament class, with its blend weight alpha.
# How these entries were fixed is documented in part1_dissertation_57/MANUAL.md, Section 5.
FINAL_CONFIGURATION = {
    "ALL": ("ET-2", 1.20),
    "PLL": ("ET-2", 1.20),
    "LF":  ("ET-1", 1.00),
    "ISL": ("ET-4", 1.20),
    "SSL": ("ET-1", 1.00),
    "ITL": ("ET-4", 1.00),
    "CL":  ("ET-2", 1.10),
}


# -------------------------------------------------------------------------- context
@dataclass
class Context:
    """Quantities learned from the training vertebrae only."""
    thresholds: dict
    bench_norm: dict
    kd_cache: dict = field(default_factory=dict)

    def zone(self, verts, lig):
        return gcz.candidate_vertices(verts, lig, self.thresholds, third_fallback=False)


def benchmark_prior(train: dict) -> dict:
    """Mean AABB-normalised expert landmark centroid per class over the training vertebrae."""
    out = {}
    for lig in LIGAMENTS:
        vals = []
        for d in train.values():
            g = (d.get("gt") or {}).get(lig)
            if g is None or d.get("verts") is None:
                continue
            v = d["verts"]
            vmin = v.min(0)
            span = v.max(0) - vmin
            span[span == 0] = 1.0
            vals.append((np.atleast_2d(g).mean(0) - vmin) / span)
        out[lig] = np.mean(vals, axis=0) if vals else np.array([0.5, 0.5, 0.5])
    return out


def make_context(train: dict) -> Context:
    return Context(thresholds=gcz.build_thresholds(train), bench_norm=benchmark_prior(train))


# ------------------------------------------------------------------------- features
def _zone_features(ctx, verts, pts, lig, lv):
    verts = np.asarray(verts, float)
    pts = np.atleast_2d(pts).astype(float)
    vmin = verts.min(0)
    vspan = verts.max(0) - vmin
    vspan[vspan == 0] = 1.0
    c = pts.mean(0)
    c_n = (c - vmin) / vspan
    cand = ctx.zone(verts, lig)
    if len(cand) < 2:
        cand = verts
    cand_c = cand.mean(0)
    cand_c_n = (cand_c - vmin) / vspan
    d_centroid = float(np.linalg.norm(c - cand_c))
    d_boundary = float(np.min(np.concatenate([c - cand.min(0), cand.max(0) - c])))
    yint = ctx.thresholds.get(lig, {}).get("y_n", {"lo": 0.0, "hi": 1.0})
    lo, hi = yint.get("lo", 0.0), yint.get("hi", 1.0)
    rel_y = float((c_n[1] - lo) / (hi - lo)) if (hi - lo) > 1e-6 else 0.5
    if len(cand) >= 3:
        ev = np.linalg.svd(cand - cand_c, full_matrices=False)[1]
        ev = ev / (ev.sum() + 1e-9)
        pca = [float(ev[0]), float(ev[1]), float(ev[2])]
    else:
        pca = [1.0, 0.0, 0.0]
    lvl = LVL_IDX.get(lv, 12) / float(len(LEVEL_ORDER) - 1)
    onehot = [1.0 if lig == L else 0.0 for L in LIGAMENTS]
    spread = cand.std(0) / vspan
    return np.asarray([c_n[0], c_n[1], c_n[2], abs(c_n[0] - 0.5), c_n[1], c_n[2],
                       d_centroid, d_boundary, rel_y, cand_c_n[0], cand_c_n[1], cand_c_n[2],
                       spread[0], spread[1], spread[2], float(len(cand)), *pca, lvl] + onehot, float)


def point_features(ctx, verts, pt, lig, lv, key=None):
    """36-D feature vector of one GCZ point (no annotations of the vertebra are used)."""
    base = _zone_features(ctx, verts, np.atleast_2d(pt), lig, lv)
    vmin = verts.min(0)
    vspan = verts.max(0) - vmin
    vspan[vspan == 0] = 1.0
    p = np.atleast_2d(pt).astype(float)[0]
    kd = ctx.kd_cache.get(key) if key is not None else None
    if kd is None:
        kd = cKDTree(verts)
        if key is not None:
            ctx.kd_cache[key] = kd
    _, idx = kd.query(p, k=min(30, len(verts)))
    nb = verts[np.atleast_1d(idx)]
    nb = nb - nb.mean(0)
    try:
        _, s, vt = np.linalg.svd(nb, full_matrices=False)
        normal = vt[2] / (np.linalg.norm(vt[2]) + 1e-9)
        curv = float(s[2] / (s.sum() + 1e-9))
    except Exception:
        normal, curv = np.array([0, 0, 1.0]), 0.0
    signed = float((p[0] - vmin[0]) / vspan[0] - 0.5)
    bp = vmin + ctx.bench_norm[lig] * vspan
    dirn = (bp - p) / vspan
    return np.concatenate([base, [normal[0], normal[1], normal[2], curv, signed],
                           [float(np.linalg.norm(p - bp)), dirn[0], dirn[1], dirn[2]]])


def residuals(gp, g):
    """Expert minus GCZ, Hungarian-paired; unpaired GCZ points use their nearest expert point."""
    gp = np.atleast_2d(gp).astype(float)
    g = np.atleast_2d(g).astype(float)
    d = np.linalg.norm(gp[:, None, :] - g[None, :, :], axis=2)
    ri, ci = linear_sum_assignment(d)
    out = np.zeros_like(gp)
    done = set()
    for a, b in zip(ri, ci):
        out[a] = g[b] - gp[a]
        done.add(a)
    if len(done) < len(gp):
        kd = cKDTree(g)
        for a in range(len(gp)):
            if a not in done:
                out[a] = g[kd.query(gp[a])[1]] - gp[a]
    return out


# Lateral ligament classes whose corrected points must stay on the side (left/right of the
# vertebra's x-centroid, the same midline used for mirroring) of the GCZ point they came from.
# Empty by default so Part 1 reproduces the dissertation; Part 2 sets {"ITL", "CL"}.
SIDE_LOCKED: set = set()


def project_to_zone(ctx, verts, pts, lig, ref=None):
    cand = ctx.zone(verts, lig)
    if len(cand) < 2:
        return pts
    pts = np.atleast_2d(pts).astype(float)
    if lig not in SIDE_LOCKED or ref is None:
        return cand[cKDTree(cand).query(pts)[1]]
    cx = float(np.asarray(verts, float)[:, 0].mean())
    ref = np.atleast_2d(ref).astype(float)
    cand_side = cand[:, 0] >= cx
    out = np.empty_like(pts)
    for i, p in enumerate(pts):
        same = cand[cand_side == (ref[i, 0] >= cx)]
        # no zone vertex on this side: the correction is not applied and the GCZ point is kept
        out[i] = same[cKDTree(same).query(p)[1]] if len(same) else ref[i]
    return out


# ---------------------------------------------------------------------- training set
def training_samples(ctx, keys, verts_of, gcz_of, gt_of):
    """Per-point samples.  meta[i] = (patient, level, ligament, is_mirrored)."""
    X, Y, meta = [], [], []
    for key in keys:
        verts = verts_of(key)
        if verts is None:
            continue
        cx = float(verts[:, 0].mean())
        gczp, gtp = gcz_of(key) or {}, gt_of(key) or {}
        for lig in LIGAMENTS:
            gp, g = gczp.get(lig), gtp.get(lig)
            if gp is None or g is None or len(gp) == 0 or len(g) == 0:
                continue
            gp = np.atleast_2d(gp).astype(float)
            res = residuals(gp, g)
            for i in range(len(gp)):
                X.append(point_features(ctx, verts, gp[i], lig, key[1], key=key))
                Y.append(res[i])
                meta.append((key[0], key[1], lig, False))
                if lig in BILATERAL:
                    pm = gp[i].copy()
                    pm[0] = 2 * cx - pm[0]
                    rm = res[i].copy()
                    rm[0] = -rm[0]
                    X.append(point_features(ctx, verts, pm, lig, key[1], key=key))
                    Y.append(rm)
                    meta.append((key[0], key[1], lig, True))
    return np.array(X), np.array(Y), meta


def correction_radius(Y, meta):
    out = {}
    for lig in LIGAMENTS:
        m = [np.linalg.norm(Y[i]) for i in range(len(meta)) if meta[i][2] == lig]
        out[lig] = float(np.mean(m) + 2 * np.std(m)) if m else 5.0
    return out


def clip(dp, radius):
    mag = float(np.linalg.norm(dp))
    return dp * (radius / mag) if mag > radius and mag > 1e-9 else dp


def blend(ctx, gp, corr, alpha, verts, lig):
    gp = np.atleast_2d(gp).astype(float)
    corr = np.atleast_2d(corr).astype(float)
    out = gp + float(alpha) * (corr - gp)
    return project_to_zone(ctx, verts, out, lig, ref=gp)


def corrected_points(ctx, model, radius, key, verts, gp, lig):
    """GCZ points of one ligament after clipped correction and zone projection (alpha = 1)."""
    gp = np.atleast_2d(gp).astype(float)
    feats = np.array([point_features(ctx, verts, pt, lig, key[1], key=key) for pt in gp])
    dps = model.predict(feats)
    corr = np.array([gp[i] + clip(np.asarray(dps[i], float), radius) for i in range(len(gp))])
    return project_to_zone(ctx, verts, corr, lig, ref=gp)


def mae(pred, gt):
    if pred is None or gt is None or len(pred) == 0 or len(gt) == 0:
        return None
    return float(gcz.hungarian_distances(pred, gt).mean())


# ---------------------------------------------------- leave-one-group-out evaluation
def held_out_records(ctx, X, Y, meta, model_name, radius, keys, verts_of, gcz_of, gt_of,
                     group_of=None, ctx_without=None, radius_without=None, original_rows_only=False):
    """Leave-one-group-out correction records used for model and alpha selection.

    group_of: maps a vertebra key to its group; default = the vertebra itself
              (leave-one-vertebra-out, dissertation); patient-level LOSO uses the subject.
    ctx_without / radius_without: optional functions returning zones/prior and correction
              radii computed without the held-out group (used in the patient-level analysis).
    original_rows_only: evaluate only the non-mirrored rows of the held-out group
              (patient-level analysis).  The dissertation run used the first len(GCZ points)
              rows of each held-out ligament, which for bilateral classes interleave
              original and mirrored rows; this is kept as published for Part 1.
    Returns {ligament: [(key, gcz_points, corrected_points, expert_points, verts), ...]}.
    """
    group_of = group_of or (lambda k: k)
    rec = {lig: [] for lig in LIGAMENTS}
    groups = list(dict.fromkeys(group_of(k) for k in keys))
    for grp in groups:
        tr = [i for i in range(len(meta)) if group_of((meta[i][0], meta[i][1])) != grp]
        if len(tr) < 20:
            continue
        model = fit(model_name, X[tr], Y[tr])
        hctx = ctx_without(grp) if ctx_without else ctx
        hrad = radius_without(grp) if radius_without else radius
        for key in [k for k in keys if group_of(k) == grp]:
            rows = [i for i in range(len(meta)) if (meta[i][0], meta[i][1]) == key]
            verts = verts_of(key)
            by_lig: dict = {}
            for i in rows:
                if not (original_rows_only and meta[i][3]):
                    by_lig.setdefault(meta[i][2], []).append(i)
            for lig, idxs in by_lig.items():
                gp, gt = (gcz_of(key) or {}).get(lig), (gt_of(key) or {}).get(lig)
                if gp is None or gt is None:
                    continue
                gp = np.atleast_2d(gp).astype(float)
                use = idxs[:len(gp)]
                feats = (np.array([point_features(hctx, verts, gp[k], lig, key[1], key=key) for k in range(len(use))])
                         if hctx is not ctx else X[use])
                dps = model.predict(feats)
                corr = gp.copy()
                for k in range(len(use)):
                    corr[k] = gp[k] + clip(np.asarray(dps[k], float), float(hrad.get(lig, np.inf)))
                rec[lig].append((key, gp, project_to_zone(hctx, verts, corr, lig, ref=gp), np.atleast_2d(gt), verts))
    return rec


def record_error(ctx, recs, alpha, lig):
    errs = [e for (_k, gp, corr, gt, verts) in recs
            if (e := mae(blend(ctx, gp, corr, alpha, verts, lig), gt)) is not None]
    return float(np.mean(errs)) if errs else float("nan")
