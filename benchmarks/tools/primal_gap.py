#!/usr/bin/env python3
"""Primal quality of two sides from run.py result JSONs (final incumbents, not times).

Per instance the primal gap of the final objective z to the reference objective z* (run.py's
"ref_objective": the MIPLIB .solu value or the reference solvers' consensus):

    gap = 0 if z == z*, 1 if there is no solution or z and z* differ in sign, else |z - z*| / max(|z|, |z*|)

(Berthold, "Measuring the impact of primal heuristics", Oper. Res. Letters 41 (2013)), averaged over
the seeds of a side. Reports the mean gap per side, the instances where B ends better / worse than A
(gap difference > --tol), and the largest changes. Files are paired by position (seed order):

    benchmarks/tools/primal_gap.py --a A-s1.json A-s2.json --b B-s1.json B-s2.json [--tol 1e-6]

(ab.py writes one run.py result per side and seed, tagged ab-A-s<seed> / ab-B-s<seed>.)
"""

import argparse
import json
from collections import defaultdict
from pathlib import Path


def pgap(z, ref):
    if z is None or ref is None:
        return 1.0
    if abs(z - ref) <= 1e-9 * max(1.0, abs(ref)):
        return 0.0
    if z * ref < 0:
        return 1.0
    return abs(z - ref) / max(abs(z), abs(ref))


def load(files):
    gaps = defaultdict(list)
    for f in files:
        for r in json.loads(Path(f).read_text())["results"]:
            if r.get("ref_objective") is None and r.get("ref_status") != "optimal":
                continue  # no reference value (infeasible or unknown)
            z = r.get("objective") if r.get("class") not in ("error", "WRONG") else None
            gaps[r["instance"]].append(pgap(z, r.get("ref_objective")))
    return {n: sum(g) / len(g) for n, g in gaps.items()}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--a", nargs="+", required=True, help="side A result JSONs (one per seed)")
    ap.add_argument("--b", nargs="+", required=True, help="side B result JSONs (one per seed)")
    ap.add_argument("--tol", type=float, default=1e-6, help="gap difference that counts as a change")
    ap.add_argument("--top", type=int, default=10)
    args = ap.parse_args()
    a, b = load(args.a), load(args.b)
    names = sorted(set(a) & set(b))
    if not names:
        raise SystemExit("no common instances")
    better = [n for n in names if b[n] < a[n] - args.tol]
    worse = [n for n in names if b[n] > a[n] + args.tol]
    ma = sum(a[n] for n in names) / len(names)
    mb = sum(b[n] for n in names) / len(names)
    print(f"instances {len(names)}  seeds A {len(args.a)} B {len(args.b)}")
    print(f"  mean primal gap  A {ma:.4f}  B {mb:.4f}")
    for side, d in (("A", a), ("B", b)):
        print(f"  {side}: at the reference {sum(d[n] <= 1e-9 for n in names)}, within 1% {sum(d[n] <= 0.01 for n in names)}, "
              f"gap 1 (no solution, or the sign differs) {sum(d[n] >= 1.0 for n in names)}")
    print(f"  B better on {len(better)}, worse on {len(worse)}, equal on {len(names) - len(better) - len(worse)}")
    ch = sorted(names, key=lambda n: b[n] - a[n])
    print("largest improvements (gap A -> B):")
    for n in [n for n in ch if b[n] < a[n] - args.tol][:args.top]:
        print(f"  {n:34s} {a[n]:.4f} -> {b[n]:.4f}")
    print("largest deteriorations:")
    for n in [n for n in ch[::-1] if b[n] > a[n] + args.tol][:args.top]:
        print(f"  {n:34s} {a[n]:.4f} -> {b[n]:.4f}")


if __name__ == "__main__":
    main()
