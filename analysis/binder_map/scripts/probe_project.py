#!/usr/bin/env python3
"""Task 9: project the full stage-3 search cost from the two-point scaling probe.

This is `probe_scaling.sh`'s own formula (its lines 31-33), applied instead of
eyeballed. Two probe arms time 2,000 queries against a 50k-target sub-DB and a
150k-target sub-DB; the ratio gives the scaling exponent, which extrapolates to
the full 500k x 500k search.

Used as the second half of the stage-2 -> stage-3 gate, so it fails closed: a
probe arm that did not run leaves an empty or absent timing, and that must stop
the chain rather than release a 512-task, ~2,100 CPU-h array.

Usage:
    python scripts/probe_project.py --t50 189 --t150 570 --ceiling 6000
"""
from __future__ import annotations

import argparse
import math

#: ~3x the plan's 2,100 CPU-h estimate. A "something is badly wrong" brake,
#: not a budget.
DEFAULT_CEILING_CPU_H = 6000.0


def valid_timing(seconds: float) -> bool:
    """A usable probe timing: finite and strictly positive.

    Zero, negative and NaN all mean the arm did not measure anything. The same
    predicate guards the CPU-h ceiling, where a non-finite value would disable
    the comparison instead of tightening it.
    """
    try:
        value = float(seconds)
    except (TypeError, ValueError):
        return False
    return math.isfinite(value) and value > 0


def project_cpu_h(t50: float, t150: float, n_full: int = 500_000,
                  n_probe_targets: int = 50_000, n_probe_queries: int = 2000,
                  threads: int = 16) -> tuple[float, float]:
    """Return (p, cpu_h) for the full search.

    p = log(t150 / t50) / log(3)   -- 3 because 150k / 50k = 3
    wall_s = t50 * (n_full / n_probe_targets) ** p * (n_full / n_probe_queries)
    cpu_h  = wall_s * threads / 3600

    Raises ValueError on a non-positive or non-finite timing rather than
    returning a number computed from garbage.
    """
    if not valid_timing(t50) or not valid_timing(t150):
        raise ValueError(f"probe timings must be finite and positive, "
                         f"got t50={t50!r} t150={t150!r}")
    p = math.log(t150 / t50) / math.log(3)
    wall_s = t50 * (n_full / n_probe_targets) ** p * (n_full / n_probe_queries)
    cpu_h = wall_s * threads / 3600
    return p, cpu_h


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--t50", type=float, required=True,
                        help="wall seconds, 2000 queries vs 50k targets")
    parser.add_argument("--t150", type=float, required=True,
                        help="wall seconds, 2000 queries vs 150k targets")
    parser.add_argument("--ceiling", type=float, default=DEFAULT_CEILING_CPU_H,
                        help=f"abort above this many CPU-h "
                             f"(default {DEFAULT_CEILING_CPU_H:g})")
    args = parser.parse_args(argv)

    # The ceiling is the brake. `--ceiling nan` (reachable through
    # GATE_CPUH_CEILING, which the sbatch forwards and `sbatch --export=ALL`
    # propagates) makes every `cpu_h > ceiling` comparison False and disables it
    # silently, so an unusable ceiling must fail the gate rather than widen it.
    if not valid_timing(args.ceiling):
        print(f"[projection] ceiling={args.ceiling!r}")
        print("  FAIL: the CPU-h ceiling is zero, negative or not a number -- "
              "it would disable the brake instead of applying it")
        return 1

    if not valid_timing(args.t50) or not valid_timing(args.t150):
        print(f"[projection] t50={args.t50!r} t150={args.t150!r}")
        print("  FAIL: a probe timing is zero, negative or not a number -- "
              "an arm did not run")
        return 1

    p, cpu_h = project_cpu_h(args.t50, args.t150)
    print(f"[projection] t50={args.t50:g}s t150={args.t150:g}s -> "
          f"p={p:.3f}, full search ~ {cpu_h:,.0f} CPU-h "
          f"(ceiling {args.ceiling:,.0f})")
    if cpu_h > args.ceiling:
        print(f"  FAIL: projected {cpu_h:,.0f} CPU-h exceeds the "
              f"{args.ceiling:,.0f} CPU-h ceiling")
        return 1
    print("  PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
