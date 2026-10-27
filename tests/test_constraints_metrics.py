"""Unit tests for the geometric constraints of GCZ-RL and the evaluation metrics (synthetic data).

Run: python -m unittest discover -s tests
"""
import sys
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code" / "ligament"))
import gcz  # noqa: E402
import gczrl  # noqa: E402
from metrics import pair_metrics  # noqa: E402


def grid_mesh(n=21):
    """Vertices of a flat 20 x 20 mm plate, x in [-10, 10], y in [0, 20], z = 0."""
    x, y = np.meshgrid(np.linspace(-10, 10, n), np.linspace(0, 20, n))
    return np.c_[x.ravel(), y.ravel(), np.zeros(n * n)]


class FixedZoneContext(gczrl.Context):
    """Context whose candidate region is a fixed vertex subset."""

    def __init__(self, zone):
        super().__init__(thresholds={}, bench_norm={})
        self._zone = zone

    def zone(self, verts, lig):
        return self._zone


class TestClip(unittest.TestCase):
    def test_long_correction_is_capped_and_keeps_direction(self):
        d = np.array([3.0, 4.0, 0.0])
        c = gczrl.clip(d, 2.5)
        self.assertAlmostEqual(np.linalg.norm(c), 2.5)
        np.testing.assert_allclose(c / np.linalg.norm(c), d / np.linalg.norm(d))

    def test_short_correction_is_unchanged(self):
        d = np.array([0.3, -0.4, 0.0])
        np.testing.assert_array_equal(gczrl.clip(d, 2.5), d)


class TestProjection(unittest.TestCase):
    def setUp(self):
        self.verts = grid_mesh()
        self.zone = self.verts[self.verts[:, 1] >= 15]   # admissible region: y >= 15
        self.ctx = FixedZoneContext(self.zone)
        self.side_locked = set(gczrl.SIDE_LOCKED)

    def tearDown(self):
        gczrl.SIDE_LOCKED = self.side_locked

    def test_output_is_a_vertex_of_the_region(self):
        pts = np.array([[1.2, 3.0, 4.0], [-7.7, 30.0, -2.0]])
        out = gczrl.project_to_zone(self.ctx, self.verts, pts, "ALL")
        zone = {tuple(v) for v in self.zone}
        self.assertTrue(all(tuple(p) in zone for p in out))

    def test_side_constraint_keeps_the_side_of_the_gcz_point(self):
        gczrl.SIDE_LOCKED = {"ITL"}
        ref = np.array([[-6.0, 16.0, 0.0]])            # GCZ point on the -x side
        target = np.array([[8.0, 16.0, 0.0]])          # correction pointing to the +x side
        free = gczrl.project_to_zone(self.ctx, self.verts, target, "ALL", ref=ref)
        locked = gczrl.project_to_zone(self.ctx, self.verts, target, "ITL", ref=ref)
        self.assertGreater(free[0, 0], 0)
        self.assertLess(locked[0, 0], 0)

    def test_side_constraint_keeps_gcz_point_when_side_is_empty(self):
        gczrl.SIDE_LOCKED = {"ITL"}
        ctx = FixedZoneContext(self.zone[self.zone[:, 0] > 0.5])   # region only on the +x side
        ref = np.array([[-6.0, 16.0, 0.0]])
        out = gczrl.project_to_zone(ctx, self.verts, np.array([[5.0, 16.0, 0.0]]), "ITL", ref=ref)
        np.testing.assert_array_equal(out, ref)

    def test_alpha_zero_returns_the_gcz_point(self):
        gp = self.zone[[3, 40]]
        corr = self.zone[[10, 50]]
        np.testing.assert_array_equal(gczrl.blend(self.ctx, gp, corr, 0.0, self.verts, "ALL"), gp)


class TestCandidateRegion(unittest.TestCase):
    def test_band_rule_and_fallback(self):
        verts = grid_mesh()
        thr = {"ALL": {"y_n": {"lo": 0.8, "hi": 1.0}, "x_half": {"lo": 0.0, "hi": 0.5},
                       "z_n": {"lo": float("nan"), "hi": float("nan")}}}
        region = gcz.candidate_vertices(verts, "ALL", thr)
        self.assertTrue(np.all(region[:, 1] >= 16.0 - 1e-9))
        thr["ALL"]["y_n"] = {"lo": 0.999, "hi": 1.0}   # fewer than 64 vertices: falls back to fewer bands
        self.assertGreaterEqual(len(gcz.candidate_vertices(verts, "ALL", thr)), 21)


class TestMetrics(unittest.TestCase):
    def test_matched_and_matching_independent_metrics(self):
        gt = np.array([[0.0, 0, 0], [10.0, 0, 0], [20.0, 0, 0], [30.0, 0, 0]])
        pred = np.array([[1.0, 0, 0], [21.0, 0, 0]])
        m = pair_metrics(pred, gt)
        self.assertEqual((m["n_pred"], m["n_ref"], m["n_matched"]), (2, 4, 2))
        self.assertAlmostEqual(m["mae_mm"], 1.0)               # Hungarian pairs 0<->1, 20<->21
        self.assertAlmostEqual(m["pass5"], 100.0)
        self.assertAlmostEqual(m["ref_coverage"], 0.5)
        self.assertAlmostEqual(m["pred_to_ref_mean_mm"], 1.0)
        self.assertAlmostEqual(m["ref_to_pred_mean_mm"], (1 + 9 + 1 + 9) / 4)
        self.assertAlmostEqual(m["ref_pass5_all"], 50.0)

    def test_missing_prediction_is_counted(self):
        m = pair_metrics(None, np.zeros((3, 3)))
        self.assertTrue(m["missing_prediction"])
        self.assertEqual(m["ref_pass5_all"], 0.0)

    def test_hungarian_is_one_to_one(self):
        d = gcz.hungarian_distances(np.array([[0.0, 0, 0], [0.1, 0, 0]]), np.array([[0.0, 0, 0], [5.0, 0, 0]]))
        self.assertEqual(len(d), 2)
        self.assertAlmostEqual(sorted(d)[1], 4.9)


if __name__ == "__main__":
    unittest.main()
