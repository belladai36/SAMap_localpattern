# LocalSAMap

LocalSAMap is an open, species-agnostic toolkit for statistical calibration and
local interpretation of cross-species cell-mapping results. It interoperates
with [SAMap](https://github.com/atarashansky/SAMap) while keeping the numerical
and inferential components usable with any compatible cross-species graph.

## What this adds

- cell-type scores aggregated from the original cross-species cell graph;
- Monte Carlo permutation p-values and Benjamini-Hochberg q-values;
- bootstrap recovery and uncertainty summaries;
- mapping specificity (margin and normalized entropy);
- module-level and homolog-level contributions;
- sensitivity analysis across ortholog-only and ortholog-plus-paralog graphs;
- optional coordinate-agreement scores for systems with defensibly comparable
  spatial or temporal coordinates.

```text
cross-species graph + annotations -> calibration and uncertainty
                                  -> local gene/module evidence
                                  -> homology-graph sensitivity
```

## Quick start

```bash
python -m pip install -e .
python -m unittest discover -s tests -v
python scripts/validate_guarantees.py
```

The numerical core is implemented in `src/localsamap/statistics.py`.

## Guarantees and limits

The exact statements, assumptions, and proofs are in
[`docs/theory.md`](docs/theory.md). In brief:

- a weighted cosine score is bounded in `[-1, 1]`;
- the plus-one Monte Carlo permutation p-value is valid under the stated
  exchangeability null;
- its simulation error has a finite-sample Hoeffding bound;
- bootstrap recovery has an analogous concentration bound;
- BH controls FDR under independence or PRDS assumptions.

These are statistical guarantees conditional on assumptions. They do not
guarantee correct cell-type homology, correct annotations, valid batch
correction, or complete homolog detection.

The deterministic simulation results are published in
[`validation/guarantee_report.md`](validation/guarantee_report.md).

## Required citations

1. Tarashansky AJ et al. (2021). *Mapping single-cell atlases throughout
   Metazoa unravels cell type evolution*. eLife 10:e66747.
   https://doi.org/10.7554/eLife.66747
2. SAMap source code: https://github.com/atarashansky/SAMap
3. Vianello SD et al. (2025). *Deconstructing the common anteroposterior
   organisation of adult bilaterian guts*. bioRxiv 2025.07.02.662275. Related
   biological motivation for distinguishing broad conservation from precise
   regional equivalence.
   https://doi.org/10.1101/2025.07.02.662275
4. Associated adult bilaterian gut code:
   https://github.com/StefanoVianello/ASICOB_GCP_AdultBilaterianGuts

## Status

Research prototype. The statistical core is tested; users must validate each
application's sampling, annotations, homolog graph, and resampling design.
