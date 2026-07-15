import numpy as np
import unittest

from localsamap import (
    benjamini_hochberg,
    coordinate_agreement,
    hoeffding_radius,
    mapping_specificity,
    permutation_pvalue,
    weighted_cosine,
)


class StatisticsTests(unittest.TestCase):
    def test_weighted_cosine_bounds_and_identity(self):
        self.assertAlmostEqual(weighted_cosine([1, 2], [1, 2], [2, 1]), 1)
        self.assertAlmostEqual(weighted_cosine([1, 0], [-1, 0]), -1)

    def test_plus_one_permutation_pvalue_never_zero(self):
        self.assertAlmostEqual(permutation_pvalue(10, [1, 2, 3]), 0.25)

    def test_bh_known_example_and_monotonicity(self):
        q = benjamini_hochberg([0.01, 0.04, 0.03, 0.2])
        self.assertTrue(np.allclose(q, [0.04, 0.0533333333, 0.0533333333, 0.2]))

    def test_specificity_and_coordinate_agreement(self):
        result = mapping_specificity([0.1, 0.8, 0.1])
        self.assertEqual(result["best_index"], 1)
        self.assertAlmostEqual(result["margin"], 0.7)
        self.assertAlmostEqual(coordinate_agreement(0.5, 0.5), 1)

    def test_hoeffding_radius_decreases_with_n(self):
        self.assertLess(hoeffding_radius(1000), hoeffding_radius(100))


if __name__ == "__main__":
    unittest.main()
