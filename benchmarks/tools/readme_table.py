#!/usr/bin/env python3
"""Markdown comparison tables (ours vs reference solvers) for README.md.

Reads one run.py result file (our times; any class other than ok/noref counts
as unsolved) and the reference times from reference/solutions.json +
miplib_solu.json, and prints a summary table and a per-instance table.

    .venv/bin/python benchmarks/tools/readme_table.py benchmarks/results/<run>.json
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

BENCH = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BENCH))

from run import load_references, sgm  # noqa: E402

SOLVED = ("optimal", "infeasible", "unbounded", "inf_or_unbd")
NAMES = {"ours": "ours", "highs": "HiGHS", "scip": "SCIP"}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("result", help="run.py result .json")
    ap.add_argument("--reference", default=str(BENCH / "reference" / "solutions.json"))
    ap.add_argument("--miplib-solu", default=str(BENCH / "reference" / "miplib_solu.json"))
    ap.add_argument("--solvers", default="ours,highs,scip")
    ap.add_argument("--no-instances", action="store_true", help="summary table only")
    args = ap.parse_args(argv)

    d = json.loads(Path(args.result).read_text())
    recs = d["results"] if isinstance(d, dict) else d
    tl = (d.get("meta", {}).get("time_limit") if isinstance(d, dict) else None) or \
        (d.get("summary", {}).get("time_limit") if isinstance(d, dict) else None) or 120.0
    refs = load_references(args.reference, args.miplib_solu)
    solvers = [s.strip() for s in args.solvers.split(",") if s.strip()]

    table = {}  # instance -> solver -> effective time (tl = unsolved), or None if not run
    for r in recs:
        n = r["instance"]
        ok = r["class"] in ("ok", "noref")
        ref = refs.get(n, {})
        if r["class"].lower() == "wrong" and r.get("status") == "optimal" and ref.get("status") == "optimal":
            # judged against an older reference value (e.g. binkar10_1's .solu, since corrected)
            z, obj = ref.get("objective"), r.get("objective")
            ok = z is not None and obj is not None and abs(obj - z) <= 1e-6 * max(1.0, abs(z))
        row = {"ours": min(r["time"], tl) if ok else tl}
        times = refs.get(n, {}).get("times", {})
        for s in solvers[1:]:
            t = times.get(s)
            if t is None or t[0] is None:
                row[s] = None
            else:
                row[s] = min(t[0], tl) if t[1] in SOLVED else tl
        table[n] = row

    print(f"| solver | solved | shifted geomean (s) |")
    print(f"|---|---|---|")
    for s in solvers:
        ts = [row[s] for row in table.values() if row.get(s) is not None]
        if not ts:
            continue
        solved = sum(1 for t in ts if t < tl)
        print(f"| {NAMES.get(s, s)} | {solved}/{len(ts)} | {sgm(ts):.1f} |")

    if args.no_instances:
        return 0
    print()
    print("| instance | " + " | ".join(NAMES.get(s, s) for s in solvers) + " |")
    print("|---" * (len(solvers) + 1) + "|")
    for n in sorted(table):
        row = table[n]
        known = [row[s] for s in solvers if row.get(s) is not None]
        best = min(known) if known else None
        cells = []
        for s in solvers:
            t = row.get(s)
            if t is None:
                cells.append("-")
            elif t >= tl:
                cells.append("TL")
            else:
                c = f"{t:.1f}"
                cells.append(f"**{c}**" if best is not None and t <= best * 1.0001 else c)
        print(f"| {n} | " + " | ".join(cells) + " |")
    return 0


if __name__ == "__main__":
    sys.exit(main())
