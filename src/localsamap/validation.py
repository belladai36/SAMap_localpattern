"""Simulation-based validation of LocalSAMap's conditional guarantees."""

from __future__ import annotations

from dataclasses import asdict, dataclass

import numpy as np

from .statistics import benjamini_hochberg, hoeffding_radius, permutation_pvalue, weighted_cosine


@dataclass(frozen=True)
class ValidationResult:
    name: str
    estimate: float
    threshold: float
    acceptance_limit: float
    passed: bool
    simulations: int
    interpretation: str

    def as_dict(self):
        return asdict(self)


def _upper_wilson(successes: int, trials: int, z: float = 2.575829) -> float:
    """One-sided upper Wilson limit; default z is conservative (99% two-sided)."""
    p = successes / trials
    z2 = z * z
    center = (p + z2 / (2 * trials)) / (1 + z2 / trials)
    radius = z * np.sqrt(p * (1 - p) / trials + z2 / (4 * trials**2)) / (1 + z2 / trials)
    return float(center + radius)


def validate_cosine_bound(rng, simulations: int = 20_000, dimension: int = 40):
    maximum = 0.0
    for _ in range(simulations):
        x, y = rng.normal(size=(2, dimension))
        w = rng.lognormal(size=dimension)
        maximum = max(maximum, abs(weighted_cosine(x, y, w)))
    return ValidationResult("weighted_cosine_bound", maximum, 1.0, 1.0, maximum <= 1.0,
                            simulations, "Maximum absolute simulated weighted cosine.")


def validate_permutation_type1(rng, simulations: int = 4_000, permutations: int = 999,
                               alpha: float = 0.05):
    """Exact exchangeable null represented by iid continuous statistics."""
    rejections = 0
    for _ in range(simulations):
        statistics = rng.normal(size=permutations + 1)
        p = permutation_pvalue(statistics[0], statistics[1:])
        rejections += p <= alpha
    upper = _upper_wilson(rejections, simulations)
    acceptance = alpha + 0.015
    return ValidationResult("permutation_type1", rejections / simulations, alpha,
                            acceptance, upper <= acceptance, simulations,
                            f"Rejection rate; conservative Wilson upper bound={upper:.4f}.")


def validate_bh_fdr(rng, simulations: int = 8_000, hypotheses: int = 100,
                    alpha: float = 0.05):
    """All-null independent Uniform(0,1) p-values; FDR=P(any rejection)."""
    false_discovery_proportions = np.empty(simulations)
    for i in range(simulations):
        q = benjamini_hochberg(rng.uniform(size=hypotheses))
        rejected = q <= alpha
        false_discovery_proportions[i] = float(np.any(rejected))
    estimate = float(false_discovery_proportions.mean())
    upper = _upper_wilson(int(false_discovery_proportions.sum()), simulations)
    acceptance = alpha + 0.012
    return ValidationResult("bh_fdr_all_null", estimate, alpha, acceptance,
                            upper <= acceptance,
                            simulations, f"Empirical FDR; conservative Wilson upper bound={upper:.4f}.")


def validate_hoeffding_coverage(rng, simulations: int = 20_000, n: int = 200,
                                probability: float = 0.3, delta: float = 0.05):
    radius = hoeffding_radius(n, delta)
    estimates = rng.binomial(n, probability, size=simulations) / n
    violations = int(np.count_nonzero(np.abs(estimates - probability) > radius))
    upper = _upper_wilson(violations, simulations)
    return ValidationResult("hoeffding_coverage", violations / simulations, delta, delta,
                            upper <= delta, simulations,
                            f"Bound-violation rate; conservative Wilson upper bound={upper:.4f}.")


def validate_bootstrap_concentration(rng, simulations: int = 20_000, replicates: int = 200,
                                     recovery_probability: float = 0.7, delta: float = 0.05):
    """Conditional iid Bernoulli recovery model used in Theorem 3."""
    radius = hoeffding_radius(replicates, delta)
    estimates = rng.binomial(replicates, recovery_probability, size=simulations) / replicates
    violations = int(np.count_nonzero(np.abs(estimates - recovery_probability) > radius))
    upper = _upper_wilson(violations, simulations)
    return ValidationResult("bootstrap_concentration", violations / simulations, delta, delta,
                            upper <= delta, simulations,
                            f"Bound-violation rate; conservative Wilson upper bound={upper:.4f}.")


def run_validation(seed: int = 20260715):
    rng = np.random.default_rng(seed)
    return [
        validate_cosine_bound(rng),
        validate_permutation_type1(rng),
        validate_bh_fdr(rng),
        validate_hoeffding_coverage(rng),
        validate_bootstrap_concentration(rng),
    ]
