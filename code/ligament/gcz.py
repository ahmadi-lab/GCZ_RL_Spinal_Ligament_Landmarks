"""Geometry-Calibrated Zoning (GCZ): the geometry-only landmark predictor.

This module is a line-by-line port of the original notebook implementation
(``fair_comparison.ipynb`` cells 4, 6, 7 and the calibration-side
``_gcz_one`` of the GCZ-RL development cells).  Random seeds, ICP settings,
thresholds and snapping rules are unchanged so that the archived predictions
can be regenerated and compared.

Main steps for one target vertebra:
1. Threshold calibration: per-ligament bands of AABB-normalised coordinates
   (mean +/- 2 SD over the training annotations).
2. Atlas selection: quick ICP from each atlas vertebra; the atlas with the
   smallest residual is retained.
3. Full ICP from the selected atlas, transfer of its expert landmarks.
4. Snapping of transferred landmarks to candidate vertices inside the
   calibrated zone, followed by farthest-point sampling to 10 points.
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import linear_sum_assignment
from scipy.spatial import cKDTree

from data_io import LIGAMENTS

N_TPL = 10
K_SIGMA = 2.0
MAX_SNAP_MM = 15.0
BILATERAL_LR = {"ITL", "LF"}
THRESHOLD_FEATURES = ("y_n", "x_half", "z_n")


# --------------------------------------------------------------------- helpers
def icp(src, tgt, n_sub=4000, max_iter=60, tol=1e-5, seed=42):
    """Rigid point-to-point ICP (src -> tgt). Returns R (3x3), t (3,)."""
    rng = np.random.default_rng(seed)
    s = np.asarray(src, float).copy()
    kd = cKDTree(tgt)
    r_acc = np.eye(3)
    t_acc = np.zeros(3)
    prev = np.inf
    for _ in range(max_iter):
        idx = rng.choice(len(s), min(n_sub, len(s)), replace=False)
        q = s[idx]
        _, nn = kd.query(q)
        tm = tgt[nn]
        mu_s = q.mean(0)
        mu_t = tm.mean(0)
        h = (q - mu_s).T @ (tm - mu_t)
        u, _, vt = np.linalg.svd(h)
        r = vt.T @ u.T
        if np.linalg.det(r) < 0:
            vt[-1] *= -1
            r = vt.T @ u.T
        t = mu_t - r @ mu_s
        s = s @ r.T + t
        r_acc = r @ r_acc
        t_acc = r @ t_acc + t
        err = float(np.linalg.norm(q @ r.T + t - tm, axis=1).mean())
        if abs(prev - err) < tol:
            break
        prev = err
    return r_acc, t_acc


def apply(r, t, pts):
    return np.asarray(pts, float) @ r.T + t


def fps(pts, n, seed=0):
    """Farthest-point sampling indices."""
    rng = np.random.default_rng(seed)
    pts = np.asarray(pts, float)
    n_all = len(pts)
    n = min(n, n_all)
    sel = [int(rng.integers(n_all))]
    dist = np.full(n_all, np.inf)
    for _ in range(n - 1):
        d = np.linalg.norm(pts - pts[sel[-1]], axis=1)
        dist = np.minimum(dist, d)
        sel.append(int(np.argmax(dist)))
    return np.array(sel)


def aabb_norm(pts, verts):
    verts = np.asarray(verts, float)
    vmin = verts.min(axis=0)
    vrng = verts.max(axis=0) - vmin
    vrng[vrng == 0] = 1.0
    return (np.asarray(pts, float) - vmin) / vrng


def threshold_features(vn):
    vn = np.asarray(vn, float)
    return {"y_n": vn[:, 1], "x_half": np.abs(vn[:, 0] - 0.5), "z_n": vn[:, 2]}


# ------------------------------------------------------------- calibration
def build_thresholds(train_data: dict) -> dict:
    """Per-ligament zone bands from training annotations only.

    ``train_data`` maps a case key to ``{"verts": V, "gt": {lig: pts}}``.
    """
    thresholds = {}
    for lig in LIGAMENTS:
        pool = {f: [] for f in THRESHOLD_FEATURES}
        for d in train_data.values():
            pts = (d.get("gt") or {}).get(lig)
            if pts is None or len(pts) == 0:
                continue
            ff = threshold_features(aabb_norm(pts, d["verts"]))
            for f in THRESHOLD_FEATURES:
                pool[f].append(ff[f])
        band = {}
        for f in THRESHOLD_FEATURES:
            if not pool[f]:
                band[f] = {"mu": float("nan"), "sigma": float("nan"), "lo": float("nan"), "hi": float("nan")}
                continue
            vals = np.concatenate(pool[f]).astype(float)
            mu = float(vals.mean())
            sd = float(vals.std(ddof=1)) if len(vals) > 1 else 0.0
            band[f] = {"mu": mu, "sigma": sd, "lo": max(0.0, mu - K_SIGMA * sd), "hi": min(1.0, mu + K_SIGMA * sd)}
        thresholds[lig] = band
    return thresholds


def candidate_vertices(verts, lig, thresholds, min_keep=64, third_fallback=True):
    """Vertices inside the calibrated zone of ``lig`` (with the original fallbacks)."""
    verts = np.asarray(verts, float)
    thr = thresholds.get(lig)
    if thr is None:
        return verts
    ff = threshold_features(aabb_norm(verts, verts))
    mask = np.ones(len(verts), dtype=bool)
    for f in THRESHOLD_FEATURES:
        lo, hi = thr[f]["lo"], thr[f]["hi"]
        if np.isnan(lo) or np.isnan(hi):
            continue
        mask &= (ff[f] >= lo) & (ff[f] <= hi)
    if mask.sum() < min_keep:
        mask = np.ones(len(verts), dtype=bool)
        for f in ("y_n", "x_half"):
            lo, hi = thr[f]["lo"], thr[f]["hi"]
            if np.isnan(lo) or np.isnan(hi):
                continue
            mask &= (ff[f] >= lo) & (ff[f] <= hi)
    if third_fallback and mask.sum() < min_keep:
        mask = np.ones(len(verts), dtype=bool)
        lo, hi = thr["y_n"]["lo"], thr["y_n"]["hi"]
        if not (np.isnan(lo) or np.isnan(hi)):
            mask &= (ff["y_n"] >= lo) & (ff["y_n"] <= hi)
    if mask.sum() == 0:
        return verts
    return verts[mask]


def build_fallback_template(train_data: dict, ref_key, icp_kwargs: dict) -> dict:
    """FPS-pooled landmark template on the reference vertebra (used only when the
    selected atlas lacks annotations for a ligament)."""
    ref_v = train_data[ref_key]["verts"]
    ref_gt = train_data[ref_key].get("gt") or {}
    aligned = {ref_key: {lig: ref_gt.get(lig) for lig in LIGAMENTS}}
    for key, d in train_data.items():
        if key == ref_key:
            continue
        r, t = icp(d["verts"], ref_v, **icp_kwargs)
        gt = d.get("gt") or {}
        aligned[key] = {lig: apply(r, t, gt[lig]) if gt.get(lig) is not None else None for lig in LIGAMENTS}
    kd = cKDTree(ref_v)
    template = {}
    for lig in LIGAMENTS:
        pts = [aligned[k].get(lig) for k in aligned if aligned[k].get(lig) is not None and len(aligned[k][lig])]
        if not pts:
            template[lig] = None
            continue
        pooled = np.vstack(pts)
        idx = fps(pooled, N_TPL, seed=0) if len(pooled) > N_TPL else np.arange(len(pooled))
        _, snn = kd.query(pooled[idx])
        template[lig] = ref_v[snn]
    return template


# ------------------------------------------------------------- prediction
def _snap(tv, query, min_ok):
    sd, nn = cKDTree(tv).query(query)
    ok = sd <= MAX_SNAP_MM
    return tv[nn[ok] if ok.sum() >= min_ok else nn]


def predict_from_alignment(v, r, t, atlas_verts, atlas_gt, thresholds, fallback,
                           third_fallback=True) -> dict:
    """Steps 3-4: warp atlas landmarks and snap them into the calibrated zones."""
    test_cx = float(v[:, 0].mean())
    crv_cx = float(atlas_verts[:, 0].mean())
    pred = {}
    for lig in LIGAMENTS:
        tpl_raw = atlas_gt.get(lig)
        if tpl_raw is None or len(tpl_raw) == 0:
            tpl_raw = fallback.get(lig)
        if tpl_raw is None or len(tpl_raw) == 0:
            pred[lig] = None
            continue
        tpl = np.asarray(tpl_raw, float)
        tv = candidate_vertices(v, lig, thresholds, third_fallback=third_fallback)
        if len(tv) == 0:
            pred[lig] = None
            continue
        if lig == "CL":
            wp = apply(r, t, tpl)
            wz = float(np.median(wp[:, 2]))
            parts = []
            for xs, zs in [("L", "S"), ("L", "I"), ("R", "S"), ("R", "I")]:
                xm = (wp[:, 0] >= test_cx) if xs == "L" else (wp[:, 0] < test_cx)
                zm = (wp[:, 2] >= wz) if zs == "S" else (wp[:, 2] < wz)
                qw = wp[xm & zm]
                if len(qw) == 0:
                    continue
                xt = (tv[:, 0] >= test_cx) if xs == "L" else (tv[:, 0] < test_cx)
                zt = (tv[:, 2] >= wz) if zs == "S" else (tv[:, 2] < wz)
                tvq = tv[xt & zt]
                if len(tvq) == 0:
                    continue
                parts.append(_snap(tvq, qw, 2))
            if not parts:
                for xs in ("L", "R"):
                    xm = (wp[:, 0] >= test_cx) if xs == "L" else (wp[:, 0] < test_cx)
                    xt = (tv[:, 0] >= test_cx) if xs == "L" else (tv[:, 0] < test_cx)
                    if xm.sum() == 0 or xt.sum() == 0:
                        continue
                    parts.append(_snap(tv[xt], wp[xm], 2))
            if not parts:
                pred[lig] = None
                continue
            cands = np.unique(np.vstack(parts), axis=0)
        elif lig in BILATERAL_LR:
            parts = []
            for ts, tvs in [(tpl[tpl[:, 0] >= crv_cx], tv[tv[:, 0] >= test_cx]),
                            (tpl[tpl[:, 0] < crv_cx], tv[tv[:, 0] < test_cx])]:
                if len(ts) == 0 or len(tvs) == 0:
                    continue
                parts.append(_snap(tvs, apply(r, t, ts), 2))
            if not parts:
                pred[lig] = None
                continue
            cands = np.unique(np.vstack(parts), axis=0)
        else:
            wp = apply(r, t, tpl)
            cands = np.unique(_snap(tv, wp, max(4, N_TPL)), axis=0)
        if len(cands) == 0:
            pred[lig] = None
        elif len(cands) <= N_TPL:
            pred[lig] = cands
        else:
            pred[lig] = cands[fps(cands, N_TPL, seed=42)]
    return pred


def select_atlas(v, atlas, quick_seeds, sample_rng):
    """Step 2: choose the atlas vertebra with the smallest quick-ICP residual."""
    kd = cKDTree(v)
    best_res, best_key = float("inf"), None
    for (key, av), seed in zip(atlas.items(), quick_seeds):
        rq, tq = icp(av, v, n_sub=1200, max_iter=40, tol=1e-5, seed=seed)
        idx = sample_rng.choice(len(av), min(1200, len(av)), replace=False)
        res = float(kd.query(apply(rq, tq, av[idx]))[0].mean())
        if res < best_res:
            best_res, best_key = res, key
    return best_key, best_res


def predict_case(v, atlas_data: dict, thresholds, fallback, *, quick_seeds, sample_seed,
                 full_seed, full_tol=1e-6, third_fallback=True):
    """Full GCZ prediction for one vertebra.

    ``atlas_data`` maps atlas keys to ``{"verts", "gt"}`` (training subjects only).
    """
    atlas = {k: d["verts"] for k, d in atlas_data.items()}
    rng = np.random.default_rng(sample_seed)
    best_key, best_res = select_atlas(v, atlas, quick_seeds, rng)
    av = atlas_data[best_key]["verts"]
    r, t = icp(av, v, n_sub=5000, max_iter=150, tol=full_tol, seed=full_seed)
    pred = predict_from_alignment(v, r, t, av, atlas_data[best_key].get("gt") or {},
                                  thresholds, fallback, third_fallback=third_fallback)
    return pred, best_key, best_res


def hungarian_distances(pred, gt):
    p = np.atleast_2d(np.asarray(pred, float))
    g = np.atleast_2d(np.asarray(gt, float))
    d = np.linalg.norm(p[:, None, :] - g[None, :, :], axis=2)
    ri, ci = linear_sum_assignment(d)
    return d[ri, ci]
