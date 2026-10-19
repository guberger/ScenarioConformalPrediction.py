import unittest

import numpy as np
from scipy.optimize import minimize
from scipy.stats import binom

from experiment_utils import select_matched_rank
from methods import (
    find_support_indices,
    fit_cascade_predictor,
    fit_threshold_predictor,
    minimum_enclosing_disk,
)


class DiskTests(unittest.TestCase):
    def test_known_disks_and_degenerate_support(self):
        acute = np.array([[0.0, 0.0], [2.0, 0.0], [1.0, np.sqrt(3)]])
        disk = minimum_enclosing_disk(acute)
        np.testing.assert_allclose(disk.center, [1, 1 / np.sqrt(3)])
        self.assertAlmostEqual(disk.radius_squared, 4 / 3)
        self.assertEqual(len(find_support_indices(acute, disk)), 3)
        obtuse = np.array([[-2.0, 0.0], [2.0, 0.0], [0.0, 0.2], [0.2, 0.1]])
        disk = minimum_enclosing_disk(obtuse)
        np.testing.assert_allclose(disk.center, [0, 0], atol=1e-14)
        self.assertAlmostEqual(disk.radius_squared, 4)
        np.testing.assert_array_equal(find_support_indices(obtuse, disk), [0, 1])
        square = np.array([[-1.0, -1.0], [-1.0, 1.0], [1.0, -1.0], [1.0, 1.0]])
        with self.assertRaisesRegex(RuntimeError, "Degenerate"):
            find_support_indices(square, minimum_enclosing_disk(square))
        np.testing.assert_allclose(minimum_enclosing_disk([[0, 0], [2, 0], [1, 0]]).center, [1, 0])
        self.assertEqual(minimum_enclosing_disk(np.ones((4, 2))).radius_squared, 0)

    def test_against_independent_convex_solver(self):
        rng = np.random.default_rng(92)
        for _ in range(12):
            z = rng.normal(size=(35, 2)) * [2, 0.7]
            disk = minimum_enclosing_disk(z)
            result = minimize(
                lambda q: q[2],
                [0.0, 0.0, 100.0],
                method="SLSQP",
                constraints={
                    "type": "ineq",
                    "fun": lambda q: q[2] - np.sum((z - q[:2]) ** 2, axis=1),
                },
                options={"ftol": 1e-10, "maxiter": 1000},
            )
            self.assertTrue(result.success, result.message)
            np.testing.assert_allclose(disk.center, result.x[:2], atol=2e-5)
            self.assertAlmostEqual(disk.radius_squared, result.x[2], places=6)
            self.assertTrue(disk.contains(z).all())

    def test_permutation_translation_scale_and_support(self):
        rng = np.random.default_rng(45)
        z = rng.normal(size=(60, 2))
        original = minimum_enclosing_disk(z)
        changed = minimum_enclosing_disk(z[rng.permutation(len(z))] * 0.003 + [3, -2])
        np.testing.assert_allclose(changed.center, original.center * 0.003 + [3, -2])
        self.assertAlmostEqual(changed.radius_squared, original.radius_squared * 0.003**2)
        indices = find_support_indices(z, original)
        reconstructed = minimum_enclosing_disk(z[indices])
        np.testing.assert_allclose(reconstructed.center, original.center, atol=1e-10)
        for i in indices:
            self.assertLess(
                minimum_enclosing_disk(np.delete(z, i, axis=0)).radius_squared,
                original.radius_squared - 1e-8,
            )

    def test_prescribed_padding_and_reoptimization(self):
        z = np.array(
            [
                [0.0, 0.2],
                [-5.0, 0.0],
                [5.0, 0.0],
                [0.0, 4.0],
                [-2.0, 0.0],
                [2.0, 0.0],
                [1.0, 1.0],
                [-1.0, 0.5],
                [0.1, -0.2],
                [0.2, 0.1],
            ]
        )
        cascade = fit_cascade_predictor(z, 1)
        self.assertEqual(cascade.stages[0]["support"], [1, 2])
        self.assertEqual(cascade.stages[0]["padding"], [0])
        # Padding uses the smallest label, NOT the next-largest residual ([0,4]).
        self.assertEqual(len(cascade.removed), 3)
        expected = minimum_enclosing_disk(z[[3, 4, 5, 6, 7, 8, 9]])
        np.testing.assert_allclose(cascade.disk.center, expected.center)
        self.assertTrue(cascade.disk.contains(z[[3, 4, 5, 6, 7, 8, 9]]).all())
        self.assertEqual(cascade.stages[-1]["removed"], [])

    def test_cascade_all_stages_and_zero_rounds(self):
        z = np.random.default_rng(9).normal(size=(70, 2))
        result = fit_cascade_predictor(z, 4)
        self.assertEqual(len(result.stages), 5)
        self.assertEqual(len(np.unique(result.removed)), 12)
        keep = np.arange(len(z))
        for stage in result.stages:
            fit = minimum_enclosing_disk(z[keep])
            np.testing.assert_allclose(stage["center"], fit.center)
            support = keep[find_support_indices(z[keep], fit)]
            np.testing.assert_array_equal(stage["support"], support)
            if stage["removed"]:
                padding = np.setdiff1d(keep, support)[: 3 - len(support)]
                np.testing.assert_array_equal(stage["padding"], padding)
            keep = np.setdiff1d(keep, stage["removed"])
        zero = fit_cascade_predictor(z, 0)
        self.assertEqual(len(zero.removed), 0)
        self.assertEqual(zero.disk.radius_squared, minimum_enclosing_disk(z).radius_squared)

    def test_sbcp_matches_direct_augmented_full_conformal(self):
        rng = np.random.default_rng(13)
        z = rng.normal(size=(39, 2))
        rank = 6
        disk, threshold = fit_threshold_predictor(z, rank)
        fitted = minimum_enclosing_disk(z)
        self.assertAlmostEqual(threshold, sorted(fitted.scores(z), reverse=True)[rank - 1])
        np.testing.assert_array_equal(disk.center, fitted.center)
        candidates = np.vstack([rng.normal(size=(35, 2)) * 2, [[30, 30], [0, 0]]])
        expected = []
        for candidate in candidates:
            augmented = np.vstack([z, candidate])
            scores = minimum_enclosing_disk(augmented).scores(augmented)
            # Account only for floating-point zero-score ties in this oracle.
            count = np.count_nonzero(scores[:-1] >= scores[-1] - 1e-10)
            expected.append(count >= rank)
        np.testing.assert_array_equal(disk.contains(candidates), expected)
        smaller, _ = fit_threshold_predictor(z, rank + 1)
        self.assertLessEqual(smaller.radius_squared, disk.radius_squared)

    def test_rank_selection_and_invalid_inputs(self):
        self.assertEqual(select_matched_rank(200, 0.1), (18, 5))
        rank, rounds = select_matched_rank(200, 0.1, delta=0.05)
        self.assertEqual(rank, 3 * (rounds + 1))
        self.assertLessEqual(binom.cdf(rank - 1, 200, 0.1), 0.05)
        self.assertGreater(binom.cdf(rank + 2, 200, 0.1), 0.05)
        self.assertEqual(select_matched_rank(200, rounds=0), (3, 0))
        for call in (
            lambda: select_matched_rank(10, 0.01),
            lambda: select_matched_rank(200, delta=0),
            lambda: fit_threshold_predictor(np.ones((10, 2)), 2),
            lambda: fit_cascade_predictor(np.ones((6, 2)), 1),
            lambda: minimum_enclosing_disk([[np.nan, 0]]),
        ):
            with self.assertRaises(ValueError):
                call()
