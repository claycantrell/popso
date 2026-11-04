#!/usr/bin/env python3
"""Where does HiGHS's presolve do more than ours? (v0.6 presolve design, PLAN.md 2026-11-03.)

Subcommands (all single-threaded, presolve only):
  diff <model.mps>      our presolved model (MILP_PRESOLVE_DUMP with original column names) vs HiGHS's
                        getPresolvedLp(), matched by column name: columns kept by one side only (type, length),
                        nonzeros on the common columns, row-length histograms, HiGHS's implied integers;
                        --lp also solves both LP relaxations with HiGHS.
  ablate <model.mps>    HiGHS presolve with one rule switched off at a time (presolve_rule_off bits):
                        sizes, time and (--lp) the LP bound of each presolved model. --rules=a,b limits the rules.
  implint <model.mps>   implied-integer detection prototype on a model (normally our presolved dump):
                        P = primal rule iterated to a fixed point (an equation whose other columns are all
                        (implied) integer with coefficients and sides integral multiples of a_j), PD = also the
                        dual rule (integral bounds and, in every row of x_j, the other columns (implied) integer
                        with integral ratios); prints the counts and whether the objective becomes integral.

    .venv/bin/python benchmarks/tools/presolve_diff.py diff benchmarks/miplib/neos-3216931-puriri.mps --lp
    .venv/bin/python benchmarks/tools/presolve_diff.py ablate benchmarks/miplib/app1-1.mps --lp \\
        --rules=none,sparsify,probing
"""

from __future__ import annotations

import argparse
import collections
import math
import os
import subprocess
import tempfile
import time

import highspy

RULES = {"none": 0, "forcing_row": 64, "forcing_col": 128, "free_subst": 256, "doubleton": 512, "dep_eq": 1024,
         "dep_free": 2048, "aggregator": 4096, "parallel": 8192, "sparsify": 16384, "probing": 32768,
         "enumeration": 65536, "dual_fixing": 131072, "col_stuffing": 262144}


def highs(path=None):
    h = highspy.Highs()
    h.setOptionValue("output_flag", False)
    h.setOptionValue("threads", 1)
    if path:
        h.readModel(path)
    return h


def lp_bound(lp) -> float:
    g = highs()
    g.passModel(lp)
    g.setOptionValue("solve_relaxation", True)
    g.run()
    return g.getInfo().objective_function_value


def columns(lp) -> dict:
    a = lp.a_matrix_
    st, ix = list(a.start_), list(a.index_)
    integ = [int(t) for t in lp.integrality_] if len(lp.integrality_) else [0] * lp.num_col_
    names = list(lp.col_names_)
    return {names[j]: dict(len=st[j + 1] - st[j], int=integ[j]) for j in range(lp.num_col_)}


def row_hist(lp) -> dict:
    c = collections.Counter(collections.Counter(list(lp.a_matrix_.index_)).values())
    out = collections.Counter()
    for k, v in c.items():
        out[min(k, 10)] += v
    return dict(sorted(out.items()))


def cmd_diff(args):
    with tempfile.TemporaryDirectory(prefix="pdiff_") as d:
        dump = os.path.join(d, "ours.mps")
        env = dict(os.environ, MILP_PRESOLVE_ONLY="1", MILP_PRESOLVE_DUMP=dump, MILP_PRESOLVE_NAMES="1")
        out = subprocess.run([args.solver, args.model, "--verbosity", "2", "--threads", "1"], env=env,
                             capture_output=True, text=True)
        for line in (out.stdout + out.stderr).splitlines():
            if line.startswith(("Presolve", "Clique")):
                print("  ours:", line)
        op = highs(dump).getLp()
    orig = columns(highs(args.model).getLp())
    hh = highs(args.model)
    hh.presolve()
    hp = hh.getPresolvedLp()
    H, O = columns(hp), columns(op)
    for lab, lp in (("HiGHS", hp), ("ours", op)):
        eq = sum(1 for i in range(lp.num_row_) if lp.row_lower_[i] == lp.row_upper_[i])
        print(f"{lab:5s}: {lp.num_row_} rows ({eq} equations), {lp.num_col_} cols, {len(lp.a_matrix_.value_)} nnz, "
              f"row lengths {row_hist(lp)}")

    def kinds(names, src):
        c = collections.Counter(("int" if src[n]["int"] else "cont", min(src[n]["len"], 5)) for n in names)
        return dict(sorted(c.items()))

    only_o, only_h = set(O) - set(H), set(H) - set(O)
    common = set(O) & set(H)
    print("cols kept by ours only:", len(only_o), kinds(only_o, O))
    print("cols kept by HiGHS only:", len(only_h), kinds(only_h, H))
    print("nnz ours - HiGHS on common cols:", sum(O[n]["len"] - H[n]["len"] for n in common))
    print("HiGHS implied integers (continuous in the original):",
          sum(1 for n in H if H[n]["int"] and not orig.get(n, {"int": 1})["int"]))
    if args.lp:
        print(f"LP bound: HiGHS {lp_bound(hp):.10g}, ours {lp_bound(op):.10g}")


def cmd_ablate(args):
    only = args.rules.split(",") if args.rules else list(RULES)
    for name in only:
        h = highs(args.model)
        h.setOptionValue("presolve_rule_off", RULES[name])
        t = time.perf_counter()
        h.presolve()
        t = time.perf_counter() - t
        p = h.getPresolvedLp()
        lpv = f" LP {lp_bound(p):.8g}" if args.lp else ""
        print(f"{name:13s} off: {p.num_row_} x {p.num_col_}, nnz {len(p.a_matrix_.value_)} ({t:.2f} s){lpv}",
              flush=True)


def cmd_implint(args):
    lp = highs(args.model).getLp()
    a = lp.a_matrix_
    st, ix, va = list(a.start_), list(a.index_), list(a.value_)
    n = lp.num_col_
    integ = [int(t) != 0 for t in lp.integrality_] if len(lp.integrality_) else [False] * n
    lo, up, rl, ru = list(lp.col_lower_), list(lp.col_upper_), list(lp.row_lower_), list(lp.row_upper_)
    rows = collections.defaultdict(list)
    for j in range(n):
        for p in range(st[j], st[j + 1]):
            rows[ix[p]].append((j, va[p]))

    def integral(v):
        return abs(v - round(v)) <= 1e-9 * max(1.0, abs(v))

    def run(dual: bool):
        impl = [False] * n
        changed = True
        while changed:
            changed = False
            for j in range(n):
                if integ[j] or impl[j] or st[j + 1] == st[j]:
                    continue
                ok_d = all(math.isinf(b) or integral(b) for b in (lo[j], up[j]))
                ok_p = False
                for p in range(st[j], st[j + 1]):
                    i, aj = ix[p], va[p]
                    others = all(c == j or ((integ[c] or impl[c] or lo[c] == up[c] and integral(lo[c])) and
                                            integral(v / aj)) for c, v in rows[i])
                    sides = all(math.isinf(s) or integral(s / aj) for s in (rl[i], ru[i]))
                    ok_p = ok_p or (rl[i] == ru[i] and others and sides)
                    ok_d = ok_d and others and sides
                if ok_p or (dual and ok_d):
                    impl[j] = changed = True
        return impl

    cont = sum(1 for j in range(n) if not integ[j])
    for lab, dual in (("P", False), ("PD", True)):
        impl = run(dual)
        objint = all(lp.col_cost_[j] == 0 or ((integ[j] or impl[j]) and integral(lp.col_cost_[j])) for j in range(n))
        print(f"{lab}: {sum(impl)} of {cont} continuous columns implied integer; objective integral: {objint}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("cmd", choices=["diff", "ablate", "implint"])
    ap.add_argument("model")
    ap.add_argument("--solver", default="build/popso")
    ap.add_argument("--lp", action="store_true", help="also solve the LP relaxations (HiGHS)")
    ap.add_argument("--rules", help="ablate: comma-separated rule names (default: all)")
    args = ap.parse_args()
    {"diff": cmd_diff, "ablate": cmd_ablate, "implint": cmd_implint}[args.cmd](args)


if __name__ == "__main__":
    main()
