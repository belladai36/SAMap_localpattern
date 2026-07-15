#!/usr/bin/env python3
"""Run deterministic simulations and write a Markdown validation report."""

from __future__ import annotations

import argparse
from pathlib import Path

from localsamap.validation import run_validation


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed", type=int, default=20260715)
    parser.add_argument("--output", type=Path, default=Path("validation/guarantee_report.md"))
    args = parser.parse_args()
    results = run_validation(args.seed)
    lines = [
        "# Simulation validation report",
        "",
        f"Random seed: `{args.seed}`",
        "",
        "| Check | Estimate | Theoretical level | Acceptance limit | Simulations | Result |",
        "|---|---:|---:|---:|---:|---|",
    ]
    for r in results:
        lines.append(f"| {r.name} | {r.estimate:.6f} | {r.threshold:.6f} | "
                     f"{r.acceptance_limit:.6f} | "
                     f"{r.simulations} | {'PASS' if r.passed else 'FAIL'} |")
    lines += ["", "## Interpretation", ""]
    lines += [f"- **{r.name}:** {r.interpretation}" for r in results]
    lines += [
        "",
        "These simulations test implementations under the assumptions stated in",
        "`docs/theory.md`. They do not establish that a biological dataset satisfies",
        "exchangeability, independence/PRDS, correct annotation, or correct homology.",
        "Permutation and FDR checks use a prespecified finite-simulation tolerance;",
        "a PASS is evidence of consistency with the theorem, not proof of it.",
        "",
    ]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text("\n".join(lines), encoding="utf-8")
    for r in results:
        print(f"{'PASS' if r.passed else 'FAIL'} {r.name}: {r.estimate:.6f}")
    raise SystemExit(0 if all(r.passed for r in results) else 1)


if __name__ == "__main__":
    main()
