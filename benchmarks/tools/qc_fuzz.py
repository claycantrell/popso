#!/usr/bin/env python3
"""Random small convex mixed-integer QCQPs and SOCPs (quadratic rows, LP format) solved by popso and SCIP; their
statuses and objectives must agree (relative 1e-5 with --gap 1e-6 given to popso by default).

    .venv/bin/python benchmarks/tools/qc_fuzz.py [--kind qcqp|soc] [--n 150] [--seed 1] [--popso build/popso]
                                                [--popso-args "--deterministic"] [--out DIR]

qcqp: 3-9 columns (continuous / integer / binary), 1-3 convex quadratic rows (sums of squares of integer linear
forms), a convex quadratic objective in a quarter of the cases. soc: second-order cones written as quadratic rows
(x'x - t^2 <= 0 with t >= 0, the same as t^2 - x'x >= 0, rotated x'x <= c u v), with shifts and constants, mixed with
linear rows and integers. Where SCIP's optimum is lower, its point is verified with lib/mps.py on popso's MPS copy of
the model: at a cone's apex SCIP's 1e-6 tolerance on the squared form lets a column reach about 1e-3, and such a point
fails the norm-form check; those cases are reported but not counted. Exit status 1 on any counted mismatch.
"""

from __future__ import annotations

import argparse
import os
import random
import subprocess
import sys
import tempfile
from pathlib import Path

BENCH = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BENCH))


def term(c, a, b=None):
    s = "+" if c >= 0 else "-"
    c = abs(c)
    if b is None:
        return f" {s} {c:.6g} {a}"
    if a == b:
        return f" {s} {c:.6g} {a}^2"
    return f" {s} {c:.6g} {a} * {b}"


def bracket(t):
    t = t.strip()
    return " " + (t[2:] if t.startswith("+ ") else t)


def gen_qcqp(rng, k):
    n = rng.randint(3, 9)
    names = [f"x{j}" for j in range(n)]
    kind = [rng.choice("CIB") for _ in range(n)]
    obj = [rng.randint(-5, 5) for _ in range(n)]
    lines = ["Minimize", " obj:" + "".join(term(c, names[j]) for j, c in enumerate(obj))]
    if k % 4 == 3:  # a convex quadratic objective too
        F = [[rng.choice([0, 0, rng.randint(-2, 2)]) for _ in range(n)] for _ in range(2)]
        S = {}
        for row in F:
            for a in range(n):
                for b in range(a, n):
                    v = row[a] * row[b] * (1 if a == b else 2)
                    if v:
                        S[(a, b)] = S.get((a, b), 0) + v
        if S:
            lines[1] += " + [" + bracket("".join(term(2 * v, names[a], names[b]) for (a, b), v in S.items())) + " ] / 2"
    lines.append("Subject To")
    for r in range(rng.randint(1, 3)):
        F = [[rng.choice([0, rng.randint(-3, 3)]) for _ in range(n)] for _ in range(rng.randint(1, 3))]
        S = {}
        for row in F:
            for a in range(n):
                for b in range(a, n):
                    v = row[a] * row[b] * (1 if a == b else 2)
                    if v:
                        S[(a, b)] = S.get((a, b), 0) + v
        if not S:
            continue
        lin = "".join(term(rng.randint(-3, 3), names[j]) for j in range(n) if rng.random() < 0.4)
        quad = " + [" + bracket("".join(term(v, names[a], names[b]) for (a, b), v in S.items())) + " ]"
        lines.append(f" q{r}: 0 x0{lin}{quad} <= {rng.randint(5, 60)}")
    for r in range(rng.randint(0, 3)):
        lines.append(f" l{r}: " + "".join(term(rng.randint(-4, 4), names[j]) for j in range(n)) + f" <= {rng.randint(0, 20)}")
    lines.append("Bounds")
    for j in range(n):
        if kind[j] == "B":
            continue
        lo, hi = (-5, 5) if kind[j] == "I" else (-10, 10)
        lines.append(f" {lo} <= {names[j]} <= {hi}")
    ints = [names[j] for j in range(n) if kind[j] == "I"]
    bins = [names[j] for j in range(n) if kind[j] == "B"]
    if ints:
        lines += ["General", " " + " ".join(ints)]
    if bins:
        lines += ["Binary", " " + " ".join(bins)]
    lines.append("End")
    return "\n".join(lines) + "\n"


def gen_soc(rng, k):
    n = rng.randint(3, 8)
    names = [f"x{j}" for j in range(n)]
    kind = [rng.choice("CCIB") for _ in range(n)]
    cones = rng.randint(1, 3)
    tnames = [f"t{i}" for i in range(cones)] + [f"v{i}" for i in range(cones)]
    obj = {nm: rng.randint(-5, 5) for nm in names}
    rows, bounds = [], {}
    for i in range(cones):
        S = rng.sample(range(n), rng.randint(1, min(4, n)))
        quad, lin, H = [], [], 0.0
        for j in S:
            c = rng.choice([1, 1, 2, 3])
            d = rng.choice([0, 0, rng.randint(-3, 3)])
            quad.append(term(c, names[j], names[j]))
            if d: lin.append(term(2 * c * d, names[j]))
            H -= c * d * d
        gamma = rng.choice([0, 0, rng.randint(1, 4)])
        H -= gamma
        form = rng.choice(["std", "std", "ge", "rot"])
        t, v = f"t{i}", f"v{i}"
        obj[t] = rng.randint(1, 4)
        bounds[t] = (0, 100)
        if form == "rot":
            obj[v] = rng.randint(1, 4)
            bounds[v] = (0, 100)
            quad.append(term(-rng.choice([1, 2, 4]), t, v))
            rows.append(f" c{i}: 0 {names[0]}{''.join(lin)} + [" + bracket("".join(quad)) + f" ] <= {H:g}")
        elif form == "ge":  # negated: c t^2 - sum ... >= -H
            ct = rng.choice([1, 2])
            neg = [term(ct, t, t)] + [q.replace(" + ", " - ", 1) if q.startswith(" + ") else q.replace(" - ", " + ", 1) for q in quad]
            nlin = [q.replace(" + ", " - ", 1) if q.startswith(" + ") else q.replace(" - ", " + ", 1) for q in lin]
            rows.append(f" c{i}: 0 {names[0]}{''.join(nlin)} + [" + bracket("".join(neg)) + f" ] >= {-H:g}")
        else:
            quad.append(term(-rng.choice([1, 2]), t, t))
            rows.append(f" c{i}: 0 {names[0]}{''.join(lin)} + [" + bracket("".join(quad)) + f" ] <= {H:g}")
    for r in range(rng.randint(1, 3)):
        rows.append(f" l{r}: " + "".join(term(rng.randint(-4, 4), nm) for nm in names) + f" >= {rng.randint(-10, 6)}")
    allv = names + [nm for nm in tnames if nm in bounds]
    lines = ["Minimize", " obj:" + "".join(term(obj.get(nm, 0), nm) for nm in allv), "Subject To"] + rows + ["Bounds"]
    for j, nm in enumerate(names):
        if kind[j] == "B":
            continue
        lo, hi = (-5, 5) if kind[j] == "I" else (-10, 10)
        lines.append(f" {lo} <= {nm} <= {hi}")
    for nm, (lo, hi) in bounds.items():
        lines.append(f" {lo} <= {nm} <= {hi}")
    ints = [names[j] for j in range(n) if kind[j] == "I"]
    bins = [names[j] for j in range(n) if kind[j] == "B"]
    if ints:
        lines += ["General", " " + " ".join(ints)]
    if bins:
        lines += ["Binary", " " + " ".join(bins)]
    lines.append("End")
    return "\n".join(lines) + "\n"


def popso_result(popso, args, path):
    p = subprocess.run([popso, "--threads", "2", "--time-limit", "60"] + args + [path], capture_output=True, text=True)
    st = [l.split(":", 1)[1].strip() for l in p.stdout.splitlines() if l.startswith("Status:")]
    ob = [float(l.split(":", 1)[1]) for l in p.stdout.splitlines()
          if l.startswith("Objective:") and l.split(":", 1)[1].strip() not in ("-", "")]
    return (st[0] if st else "ERR " + (p.stderr.strip().splitlines() or ["?"])[-1]), (ob[0] if ob else None)


def scip_point_fails_norm_form(popso, path, m):
    """True if SCIP's best point fails lib/mps.py verification (norm form for cone rows) on popso's MPS copy."""
    from lib.mps import read_mps, verify

    mps, amap = path + ".mps", path + ".map"
    subprocess.run([popso, path, "--write-anon", mps, "--anon-keep-order", "--anon-map", amap], capture_output=True)
    ren = dict(reversed(l.split()) for l in open(amap) if len(l.split()) == 2)
    s = m.getBestSol()
    vals = {ren[v.name]: m.getSolVal(s, v) for v in m.getVars() if v.name in ren}
    return not verify(read_mps(mps), vals).feasible


def main(argv=None):
    import pyscipopt

    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--kind", choices=["qcqp", "soc"], default="qcqp")
    ap.add_argument("--n", type=int, default=150)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--popso", default=str(BENCH.parent / "build" / "popso"))
    ap.add_argument("--popso-args", default="--gap 1e-6")
    ap.add_argument("--out", help="directory for the models (default: a temporary one)")
    args = ap.parse_args(argv)
    out = args.out or tempfile.mkdtemp(prefix="qc_fuzz_")
    os.makedirs(out, exist_ok=True)
    rng = random.Random(args.seed)
    gen = gen_qcqp if args.kind == "qcqp" else gen_soc
    bad = apex = 0
    for k in range(args.n):
        path = os.path.join(out, f"{args.kind}{k}.lp")
        open(path, "w").write(gen(rng, k))
        m = pyscipopt.Model()
        m.hideOutput()
        m.readProblem(path)
        m.setParam("limits/time", 60)
        m.optimize()
        st = m.getStatus()
        ref = m.getObjVal() if st == "optimal" else None
        pst, pob = popso_result(args.popso, args.popso_args.split(), path)
        ok = (st == "infeasible" and pst == "INFEASIBLE") or st == "timelimit" or (
            ref is not None and pst == "OPTIMAL" and pob is not None and abs(pob - ref) <= 1e-5 * max(1, abs(ref)))
        if ok:
            continue
        if ref is not None and pob is not None and pob > ref and scip_point_fails_norm_form(args.popso, path, m):
            apex += 1
            print(f"case {k}: SCIP {ref} below popso {pob}: SCIP's point fails the norm form (apex tolerance)")
            continue
        bad += 1
        print(f"case {k}: SCIP {st} {ref}  popso {pst} {pob}  ({path})", flush=True)
    print(f"{args.n} {args.kind} cases: {bad} mismatches, {apex} SCIP apex-tolerance cases")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
