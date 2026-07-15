import unittest

import numpy as np

from localsamap.validation import (
    validate_bh_fdr,
    validate_bootstrap_concentration,
    validate_cosine_bound,
    validate_hoeffding_coverage,
    validate_permutation_type1,
)


class ValidationSmokeTests(unittest.TestCase):
    def setUp(self):
        self.rng = np.random.default_rng(12345)

    def test_simulation_validators_execute(self):
        results = [
            validate_cosine_bound(self.rng, simulations=50),
            validate_permutation_type1(self.rng, simulations=200, permutations=99),
            validate_bh_fdr(self.rng, simulations=200, hypotheses=20),
            validate_hoeffding_coverage(self.rng, simulations=200, n=100),
            validate_bootstrap_concentration(self.rng, simulations=200, replicates=100),
        ]
        self.assertEqual(len(results), 5)
        self.assertTrue(all(0 <= r.estimate <= 1 for r in results))


if __name__ == "__main__":
    unittest.main()

