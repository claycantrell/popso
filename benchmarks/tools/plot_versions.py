#!/usr/bin/env python3
"""SVG plot "instances solved vs. time" for several popso versions on the same set.

One step line per version (oldest lightest, newest darkest), from run.py result files of
the same instance set and settings.

    .venv/bin/python benchmarks/tools/plot_versions.py \
        v0.1.0=benchmarks/results/A.json v0.2.0=benchmarks/results/B.json ... \
        --title "popso over time — MIPLIB 2017 subset" --out docs/versions.svg
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

BENCH = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BENCH))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import plot_profile as pp  # noqa: E402
from run import load_references, sgm  # noqa: E402

# ordinal blue ramp (validated): light mode light->dark, dark mode dark->light (newest most prominent)
LIGHT = ["#86b6ef", "#3987e5", "#1c5cab", "#0d366b"]
DARK = ["#184f95", "#2a78d6", "#6da7ec", "#b7d3f6"]


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("runs", nargs="+", help="LABEL=result.json, oldest first (at most 4)")
    ap.add_argument("--title", required=True)
    ap.add_argument("--subtitle", default="")
    ap.add_argument("--time-limit", type=float, default=120.0)
    ap.add_argument("--reference", default=str(BENCH / "reference" / "solutions.json"))
    ap.add_argument("--miplib-solu", default=str(BENCH / "reference" / "miplib_solu.json"))
    ap.add_argument("--out", required=True)
    args = ap.parse_args(argv)
    if len(args.runs) > len(LIGHT):
        ap.error(f"at most {len(LIGHT)} versions")
    refs = load_references(args.reference, args.miplib_solu)

    # pick ramp steps spread over the ramp, newest = darkest (light) / lightest (dark)
    k = len(args.runs)
    idx = [round(i * (len(LIGHT) - 1) / max(1, k - 1)) for i in range(k)]
    times, n = {}, None
    for i, spec in enumerate(args.runs):
        label, path = spec.split("=", 1)
        key = f"v{9 - i}"  # ties in the end labels: newest listed first
        t, m = pp.solve_times(path, ["ours"], refs)
        if n is not None and m != n:
            ap.error(f"{label}: {m} instances, expected {n} (all runs must use the same set)")
        n = m
        times[key] = t["ours"]
        g = sgm([min(x, args.time_limit) if x is not None else args.time_limit for x in t["ours"]])
        pp.NAMES[key] = f"{label} · {g:.1f} s"
        pp.COLORS[key] = (LIGHT[idx[i]], DARK[idx[i]])
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text(pp.svg(times, n, args.title, args.subtitle, args.time_limit, right=180))
    print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
