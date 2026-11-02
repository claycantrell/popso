#!/usr/bin/env python3
"""A/B comparison of two solver builds over several random seeds.

Runs benchmarks/run.py for build A and build B (interleaved per seed, so slow
drifts in machine load affect both), then reports:
  * wrong answers per side (any wrong answer fails the comparison),
  * solved counts,
  * shifted geometric mean time (shift 1 s, unsolved = time limit) per side,
  * the speed ratio B/A with a 90% bootstrap confidence interval over instances,
  * per-family ratios and the instances with the largest changes.

Usage:
  .venv/bin/python benchmarks/ab.py --a build-base/popso --b build/popso --seeds 1,2,3 --sizes small,medium -t 60
  .venv/bin/python benchmarks/ab.py --a build/popso --b build/popso --b-args "--deterministic" ...

Extra arguments for each side go in --a-args / --b-args. Instance selection
options (--sizes, --families, --max-size, positional globs, --miplib, --set) are
forwarded to run.py.
"""

import argparse
import json
import math
import random
import re
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

BENCH = Path(__file__).resolve().parent
RUN = BENCH / "run.py"


def run_side(solver, extra, seed, args, label, env_extra=""):
    cmd = [sys.executable, str(RUN), "--solver", f"{solver} --seed {seed} {extra}".strip(),
           "-t", str(args.time_limit), "--threads", str(args.threads), "-j", "1",
           "--mem-limit", str(args.mem_limit), "--gap", str(args.gap), "--tag", f"ab-{label}-s{seed}"]
    for opt in ("sizes", "families", "max_size"):
        v = getattr(args, opt)
        if v:
            cmd += ["--" + opt.replace("_", "-"), v]
    if args.miplib:
        cmd.append("--miplib")
    for st in args.set:
        cmd += ["--set", st]
    cmd += args.instances
    print(f"[ab] {label} seed {seed}: {' '.join(cmd[2:5])} ...", flush=True)
    import os
    env = dict(os.environ)
    for kv in env_extra.split():
        k, _, v = kv.partition("=")
        env[k] = v
    out = subprocess.run(cmd, capture_output=True, text=True, env=env)
    m = re.search(r"^results: (.+\.json)$", out.stdout, re.M)
    if not m:
        sys.stderr.write(out.stdout[-3000:] + out.stderr[-3000:])
        raise SystemExit(f"[ab] run.py failed for {label} seed {seed}")
    return json.loads(Path(m.group(1)).read_text())


def family(name):
    return name.split("_")[0] if "_" in name else name


def sgm(times, shift=1.0):
    return math.exp(sum(math.log(t + shift) for t in times) / len(times)) - shift


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--a", required=True, help="solver binary A (baseline)")
    ap.add_argument("--b", required=True, help="solver binary B (candidate)")
    ap.add_argument("--a-args", default="")
    ap.add_argument("--b-args", default="")
    ap.add_argument("--a-env", default="", help="environment for A, e.g. 'MILP_NO_RCFIX=1'")
    ap.add_argument("--b-env", default="")
    ap.add_argument("--seeds", default="1,2,3")
    ap.add_argument("-t", "--time-limit", type=float, default=60)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--mem-limit", type=int, default=1500)
    ap.add_argument("--gap", type=float, default=1e-6)
    ap.add_argument("--sizes")
    ap.add_argument("--families")
    ap.add_argument("--max-size")
    ap.add_argument("--miplib", action="store_true")
    ap.add_argument("--set", action="append", default=[], help="instance list (see run.py --set)")
    ap.add_argument("--out", help="write the comparison as JSON here")
    ap.add_argument("instances", nargs="*")
    args = ap.parse_args()
    seeds = [int(s) for s in args.seeds.split(",")]

    # times[side][instance] -> list over seeds
    times = {"A": defaultdict(list), "B": defaultdict(list)}
    solved = {"A": 0, "B": 0}
    wrong = {"A": [], "B": []}
    runs = 0
    for seed in seeds:
        for side, solver, extra, env in (("A", args.a, args.a_args, args.a_env), ("B", args.b, args.b_args, args.b_env)):
            res = run_side(solver, extra, seed, args, side, env)
            for r in res["results"]:
                name = r["instance"]
                ok = r["class"] == "ok"
                if r["class"] == "WRONG":
                    wrong[side].append((name, seed, r.get("message", "")))
                solved[side] += ok
                times[side][name].append(min(r["time"], args.time_limit) if ok else args.time_limit)
        runs += 1

    names = sorted(set(times["A"]) & set(times["B"]))
    if not names:
        raise SystemExit("[ab] no common instances")
    # Per-instance geometric mean over seeds, then shifted geomean over instances.
    per = {n: (sgm(times["A"][n]), sgm(times["B"][n])) for n in names}
    a_all = [per[n][0] for n in names]
    b_all = [per[n][1] for n in names]
    ga, gb = sgm(a_all), sgm(b_all)
    ratio = (gb + 1) / (ga + 1)  # ratio of shifted geomeans

    rng = random.Random(0)
    boots = []
    for _ in range(2000):
        sample = [rng.choice(names) for _ in names]
        ba = sgm([per[n][0] for n in sample])
        bb = sgm([per[n][1] for n in sample])
        boots.append((bb + 1) / (ba + 1))
    boots.sort()
    lo, hi = boots[int(0.05 * len(boots))], boots[int(0.95 * len(boots))]

    print()
    print(f"instances {len(names)} x seeds {len(seeds)}   (time limit {args.time_limit:g}s, threads {args.threads})")
    for side in ("A", "B"):
        print(f"  {side}: solved {solved[side]}/{len(names) * len(seeds)}  wrong {len(wrong[side])}  "
              f"shifted geomean {sgm(a_all if side == 'A' else b_all):.3f}s")
    verdict = "FASTER" if hi < 1.0 else "SLOWER" if lo > 1.0 else "no significant difference"
    print(f"  ratio B/A of shifted geomeans: {ratio:.3f}   90% CI [{lo:.3f}, {hi:.3f}]  -> B is {verdict}")

    fam = defaultdict(list)
    for n in names:
        fam[family(n)].append(n)
    print("\nper family (B/A ratio of shifted geomeans):")
    for f in sorted(fam):
        fa = sgm([per[n][0] for n in fam[f]])
        fb = sgm([per[n][1] for n in fam[f]])
        print(f"  {f:14s} {len(fam[f]):3d} inst  A {fa:8.3f}s  B {fb:8.3f}s  ratio {(fb + 1) / (fa + 1):.3f}")

    changes = sorted(names, key=lambda n: (per[n][1] + 1) / (per[n][0] + 1))
    print("\nlargest speedups:")
    for n in changes[:5]:
        print(f"  {n:34s} A {per[n][0]:8.2f}s  B {per[n][1]:8.2f}s")
    print("largest slowdowns:")
    for n in changes[::-1][:5]:
        print(f"  {n:34s} A {per[n][0]:8.2f}s  B {per[n][1]:8.2f}s")

    for side in ("A", "B"):
        for name, seed, msg in wrong[side]:
            print(f"!!! WRONG {side} {name} seed {seed}: {msg}")
    if args.out:
        Path(args.out).write_text(json.dumps({
            "a": args.a, "b": args.b, "a_args": args.a_args, "b_args": args.b_args, "seeds": seeds,
            "ratio": ratio, "ci90": [lo, hi], "per_instance": per, "wrong": wrong, "solved": solved}, indent=1))
    sys.exit(1 if wrong["B"] else 0)


if __name__ == "__main__":
    main()
