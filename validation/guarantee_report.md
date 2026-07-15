# Simulation validation report

Random seed: `20260715`

| Check | Estimate | Theoretical level | Acceptance limit | Simulations | Result |
|---|---:|---:|---:|---:|---|
| weighted_cosine_bound | 0.841071 | 1.000000 | 1.000000 | 20000 | PASS |
| permutation_type1 | 0.049250 | 0.050000 | 0.065000 | 4000 | PASS |
| bh_fdr_all_null | 0.052000 | 0.050000 | 0.062000 | 8000 | PASS |
| hoeffding_coverage | 0.002350 | 0.050000 | 0.050000 | 20000 | PASS |
| bootstrap_concentration | 0.002000 | 0.050000 | 0.050000 | 20000 | PASS |

## Interpretation

- **weighted_cosine_bound:** Maximum absolute simulated weighted cosine.
- **permutation_type1:** Rejection rate; conservative Wilson upper bound=0.0588.
- **bh_fdr_all_null:** Empirical FDR; conservative Wilson upper bound=0.0588.
- **hoeffding_coverage:** Bound-violation rate; conservative Wilson upper bound=0.0034.
- **bootstrap_concentration:** Bound-violation rate; conservative Wilson upper bound=0.0030.

These simulations test implementations under the assumptions stated in
`docs/theory.md`. They do not establish that a biological dataset satisfies
exchangeability, independence/PRDS, correct annotation, or correct homology.
Permutation and FDR checks use a prespecified finite-simulation tolerance;
a PASS is evidence of consistency with the theorem, not proof of it.
