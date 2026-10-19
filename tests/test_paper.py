import tempfile
import unittest

import numpy as np
import pandas as pd

from methods import find_support_indices, minimum_enclosing_disk
from paper_experiments import main, schedules


class PaperTests(unittest.TestCase):
    def test_support_detects_center_change_near_right_triangle(self):
        points = np.array([[-1.0, 0.0], [1.0, 0.0], [0.0, 1.00001]])
        support = find_support_indices(points, minimum_enclosing_disk(points))
        np.testing.assert_array_equal(support, [0, 1, 2])

    def test_exact_paper_design_and_outputs(self):
        self.assertEqual(schedules(18), [(3, 3, 3, 3, 3, 3), (9, 9)])
        with tempfile.TemporaryDirectory() as root:
            output = main(["--runs", "2", "--n-test", "60", "--horizons", "30", "--outdir", root])
            rows = pd.read_csv(output / "runs.csv")
            self.assertEqual(set(rows.method), {"hat", "check", "bar_L6", "bar_L2"})
            self.assertTrue(rows.valid.all())
            self.assertTrue((rows[rows.method == "hat"].fits == 6).all())
            self.assertTrue((rows[rows.method == "hat"].removed == 15).all())
            for name in (
                "sets.png",
                "coverage.png",
                "volume.png",
                "runtime.png",
                "report.md",
                "paired.csv",
            ):
                self.assertTrue((output / name).is_file(), name)
            with np.load(output / "first_run_T30.npz") as data:
                self.assertEqual(data["bar_L2_centers"].shape, (2, 2))
                self.assertEqual(data["bar_L6_centers"].shape, (6, 2))


if __name__ == "__main__":
    unittest.main()
