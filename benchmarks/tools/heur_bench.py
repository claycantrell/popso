#!/usr/bin/env python3
"""Standalone primal-heuristic comparison (MILP_FPR_BENCH): runs the solver binary with the
bench knob on each instance of a set (one process at a time, single-threaded, memory-capped) and
tabulates, per heuristic, whether a feasible point was found and how fast.

    .venv/bin/python benchmarks/tools/heur_bench.py --solver build/popso --set miplib_wide_core \
        --effort 1000000000 --out results.json

The solver prints one "HEURBENCH <name> found=.. time_first=.. obj_first=.. ..." line per
heuristic (see src/milp/heuristics/fpr.cpp, heuristic_bench). Times are measured inside the solver
on the presolved problem (presolve time excluded).
"""

import argparse
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

BENCH = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BENCH))
from lib.proc import run_capped  # noqa: E402


def read_set(name):
    p = Path(name)
    if not p.exists():
        p = BENCH / "sets" / f"{name}.txt"
    out = []
    for line in p.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        f = Path(line)
        if not f.is_absolute():
            f = BENCH.parent / f
        out.append(f)
    return out


def parse(text):
    res = {}
    for line in text.splitlines():
        if not line.startswith("HEURBENCH "):
            continue
        parts = line.split()
        name = parts[1]
        kv = {}
        for tok in parts[2:]:
            if "=" in tok:
                k, v = tok.split("=", 1)
                try:
                    kv[k] = float(v) if k not in ("strategy", "repairs", "leaflps") else v
                except ValueError:
                    kv[k] = v
            else:
                kv[tok] = True
        res[name] = kv
    return res


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--solver", required=True)
    ap.add_argument("--set", action="append", default=[], help="set name in benchmarks/sets or a file")
    ap.add_argument("instances", nargs="*", help="instance files (in addition to --set)")
    ap.add_argument("--effort", type=int, default=1000000000)
    ap.add_argument("--timeout", type=float, default=300.0, help="wall-clock cap per instance (s)")
    ap.add_argument("--mem-limit", type=int, default=2048, help="RSS cap per process (MB)")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--each", action="store_true", help="also run each FPR strategy alone")
    ap.add_argument("--lp", action="store_true", help="also run the LP-guided FPR after the root LP (fpr-lp)")
    ap.add_argument("--fj-from", help="skip FJ; take its results from this earlier --out JSON")
    ap.add_argument("--env", action="append", default=[], help="extra KEY=VALUE for the solver (knobs)")
    ap.add_argument("--sample", type=int, default=0, help="take every k-th instance of the list")
    ap.add_argument("--out", help="write the results as JSON")
    ap.add_argument("--resume", action="store_true", help="skip instances already in --out")
    args = ap.parse_args(argv)

    files = [Path(f) for f in args.instances]
    for s in args.set:
        files += read_set(s)
    if args.sample > 1:
        files = files[:: args.sample]
    results = {}
    if args.resume and args.out and Path(args.out).exists():
        results = json.loads(Path(args.out).read_text())
    os.environ["MILP_FPR_BENCH"] = str(args.effort)  # inherited by the solver processes
    if args.each:
        os.environ["MILP_FPR_BENCH_EACH"] = "1"
    if args.lp:
        os.environ["MILP_FPR_BENCH_LP"] = "1"
    fj_from = json.loads(Path(args.fj_from).read_text()) if args.fj_from else None
    if fj_from is not None:
        os.environ["MILP_FPR_BENCH_NOFJ"] = "1"
    for kv in args.env:
        k, v = kv.split("=", 1)
        os.environ[k] = v
    for f in files:
        name = f.name.split(".")[0]
        if name in results:
            continue
        with tempfile.TemporaryFile("w+") as out:
            cmd = [args.solver, str(f), "--threads", "1", "--seed", str(args.seed), "--verbosity", "1"]
            rc, reason, wall, peak = run_capped(cmd, args.timeout, args.mem_limit, stdout=out,
                                                stderr=subprocess.STDOUT)
            out.seek(0)
            text = out.read()
        r = parse(text)
        if fj_from is not None and "fj" in fj_from.get(name, {}):
            r["fj"] = fj_from[name]["fj"]
        r["_wall"] = wall
        r["_status"] = reason or ("ok" if rc == 0 else f"rc={rc}")
        results[name] = r
        fj, fpr = r.get("fj", {}), r.get("fpr", {})

        def show(h):
            if not h:
                return "      -     "
            return f"{'Y' if h.get('found') else '.'} {h.get('time_first', -1):8.3f}s" if h.get("found") else \
                f". ({h.get('time', 0):6.2f}s)"
        print(f"{name:32s} fj {show(fj)}   fpr {show(fpr)} {fpr.get('strategy', '')}  [{r['_status']}, "
              f"{wall:.1f}s]", flush=True)
        if args.out:
            Path(args.out).write_text(json.dumps(results, indent=1))
    names = sorted({k for r in results.values() for k in r if not k.startswith("_") and k not in ("problem", "lp")})
    print()
    for h in names:
        found = [n for n, r in results.items() if r.get(h, {}).get("found")]
        times = sorted(results[n][h]["time_first"] for n in found)
        med = times[len(times) // 2] if times else float("nan")
        print(f"{h:18s} found {len(found):3d}/{len(results)}  median time-to-first {med:.3f}s")
    fj_set = {n for n, r in results.items() if r.get("fj", {}).get("found")}
    fpr_set = {n for n, r in results.items() if r.get("fpr", {}).get("found")}
    print(f"FPR only: {sorted(fpr_set - fj_set)}")
    print(f"FJ only:  {sorted(fj_set - fpr_set)}")
    print(f"either: {len(fj_set | fpr_set)}, neither: {len(results) - len(fj_set | fpr_set)}")


if __name__ == "__main__":
    sys.exit(main())
