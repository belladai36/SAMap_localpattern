"""Small, auditable statistical primitives used by LocalSAMap."""

from __future__ import annotations

import numpy as np


def weighted_cosine(x, y, weights=None) -> float:
    """Weighted cosine similarity; weights must be finite and nonnegative."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if x.shape != y.shape:
        raise ValueError("x and y must have the same shape")
    w = np.ones_like(x) if weights is None else np.asarray(weights, dtype=float)
    if w.shape != x.shape or not np.all(np.isfinite(w)) or np.any(w < 0):
        raise ValueError("weights must match x and be finite and nonnegative")
    if not np.all(np.isfinite(x)) or not np.all(np.isfinite(y)):
        raise ValueError("x and y must be finite")
    denominator = np.sqrt(np.sum(w * x * x) * np.sum(w * y * y))
    if denominator == 0:
        raise ValueError("weighted norms must be nonzero")
    return float(np.clip(np.sum(w * x * y) / denominator, -1.0, 1.0))


def permutation_pvalue(observed: float, null_statistics, alternative="greater") -> float:
    """Plus-one Monte Carlo p-value from valid exchangeable permutations."""
    null = np.asarray(null_statistics, dtype=float)
    if null.ndim != 1 or null.size == 0 or not np.all(np.isfinite(null)):
        raise ValueError("null_statistics must be a nonempty finite vector")
    if alternative == "greater":
        exceedances = np.count_nonzero(null >= observed)
    elif alternative == "less":
        exceedances = np.count_nonzero(null <= observed)
    elif alternative == "two-sided":
        exceedances = np.count_nonzero(np.abs(null) >= abs(observed))
    else:
        raise ValueError("alternative must be greater, less, or two-sided")
    return float((exceedances + 1) / (null.size + 1))


def benjamini_hochberg(pvalues):
    """Benjamini-Hochberg adjusted p-values, preserving input shape."""
    p = np.asarray(pvalues, dtype=float)
    if np.any(~np.isfinite(p)) or np.any((p < 0) | (p > 1)):
        raise ValueError("p-values must be finite and in [0, 1]")
    flat = p.ravel()
    order = np.argsort(flat)
    ranked = flat[order]
    adjusted = ranked * flat.size / np.arange(1, flat.size + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    output = np.empty_like(adjusted)
    output[order] = np.minimum(adjusted, 1.0)
    return output.reshape(p.shape)


def hoeffding_radius(n: int, delta: float = 0.05) -> float:
    """Two-sided Hoeffding error radius at confidence 1-delta."""
    if n <= 0 or not 0 < delta < 1:
        raise ValueError("n must be positive and delta must lie in (0, 1)")
    return float(np.sqrt(np.log(2 / delta) / (2 * n)))


def mapping_specificity(scores):
    """Return best index, best-vs-second margin, and normalized entropy."""
    s = np.asarray(scores, dtype=float)
    if s.ndim != 1 or s.size < 2 or np.any(~np.isfinite(s)) or np.any(s < 0):
        raise ValueError("scores must be a finite nonnegative vector of length >= 2")
    total = s.sum()
    if total == 0:
        raise ValueError("at least one score must be positive")
    order = np.argsort(s)
    q = s / total
    positive = q > 0
    entropy = -np.sum(q[positive] * np.log(q[positive])) / np.log(s.size)
    return {"best_index": int(order[-1]), "margin": float(s[order[-1]] - s[order[-2]]),
            "normalized_entropy": float(entropy)}


def coordinate_agreement(position_a: float, position_b: float, sigma: float = 0.2) -> float:
    """Gaussian agreement for defensibly comparable coordinates in [0, 1]."""
    if not (0 <= position_a <= 1 and 0 <= position_b <= 1) or sigma <= 0:
        raise ValueError("positions must be in [0, 1] and sigma must be positive")
    return float(np.exp(-((position_a - position_b) ** 2) / (2 * sigma**2)))
