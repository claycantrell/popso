#!/usr/bin/env python3
"""Random small convex QPs and MIQPs (quadratic objective, linear constraints) solved by
popso and a reference: HiGHS for continuous QPs, SCIP for MIQPs (and for QPs where HiGHS and
popso disagree); an "unbounded" answer of popso that the references do not share is checked with
an LP certificate (feasible point + recession direction with Q d = 0, c'd < 0). Statuses and objectives must agree (relative 1e-6) and popso's solutions must pass
lib/mps.py verification (which recomputes the quadratic objective).

    .venv/bin/python benchmarks/tools/qp_fuzz.py [--n 200] [--seed 1] [--build build]
                                                [--format mps|lp] [--mip-frac 0.5] [--popso-args "..."]

Models: Q = G G^T with integer G of random rank (columns without curvature, singular Q),
free / one-sided / boxed / fixed columns, equality, ranged and one-sided rows (feasible at a
random point, sometimes made infeasible), occasional unbounded directions, maximization with
a concave objective. Runs one solver process at a time; exit status 1 on any mismatch (the
offending model is kept in --out, default: a temporary directory).
"""

from __future__ import annotations

import argparse
import re
import shlex
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

BENCH = Path(__file__).resolve().parents[1]
ROOT = BENCH.parent
sys.path.insert(0, str(BENCH))

from lib import mps, refsolvers, solfile  # noqa: E402
from lib.modeler import Model, lin  # noqa: E402

INF = float("inf")


def random_model(rng, k, mip_frac):
    m = Model(f"qpfuzz{k}")
    n = int(rng.integers(3, 13))
    mip = rng.random() < mip_frac
    x0 = rng.integers(-3, 4, n).astype(float)
    xs = []
    for j in range(n):
        kind = rng.choice(["box", "lower", "upper", "free", "fixed"], p=[0.45, 0.2, 0.1, 0.2, 0.05])
        integer = mip and rng.random() < 0.6
        lo, up = -INF, INF
        if kind == "box":
            lo, up = x0[j] - float(rng.integers(0, 4)), x0[j] + float(rng.integers(1, 4))
        elif kind == "lower":
            lo = x0[j] - float(rng.integers(0, 3))
        elif kind == "upper":
            up = x0[j] + float(rng.integers(0, 3))
        elif kind == "fixed":
            lo = up = x0[j]
        if integer and (lo == -INF or up == INF):
            lo, up = (x0[j] - 3 if lo == -INF else lo), (x0[j] + 3 if up == INF else up)
        xs.append(m.add_var(lb=lo, ub=up, vtype="I" if integer else "C", obj=float(rng.integers(-5, 6)),
                            name=f"x{j}"))
    # Curvature: Q = G G^T (rank r); columns outside the support of G are linear.
    r = int(rng.integers(1, n + 1))
    support = [j for j in range(n) if rng.random() < 0.8] or [0]
    G = np.zeros((n, r))
    for j in support:
        G[j] = rng.integers(-2, 3, r)
    Q = G @ G.T
    # Unbounded directions are possible (free linear columns); keep most models bounded by giving
    # free columns curvature.
    for j in range(n):
        if xs[j].lb == -INF and xs[j].ub == INF and Q[j, j] == 0 and rng.random() < 0.8:
            Q[j, j] = 1.0
    for i in range(int(rng.integers(0, 7))):
        idx = [j for j in range(n) if rng.random() < 0.5] or [int(rng.integers(n))]
        a = rng.integers(-4, 5, len(idx)).astype(float)
        a[a == 0] = 1.0
        act = float(a @ x0[idx])
        e = lin(a, [xs[j] for j in idx])
        t = rng.random()
        slack = float(rng.integers(0, 4))
        if rng.random() < 0.05:
            slack = -float(rng.integers(1, 4))  # infeasible at x0 (may be infeasible)
        if t < 0.2:
            m.add_constr(e == act, f"r{i}")
        elif t < 0.6:
            m.add_constr(e <= act + slack, f"r{i}")
        else:
            m.add_constr(e >= act - slack, f"r{i}")
    sign = 1.0
    if rng.random() < 0.3:
        m.sense = "max"
        sign = -1.0
        for v in xs:
            v.obj = -v.obj
    for i in range(n):
        if Q[i, i] != 0:
            m.add_quad(xs[i], xs[i], sign * 0.5 * Q[i, i])
        for j in range(i):
            if Q[i, j] != 0:
                m.add_quad(xs[i], xs[j], sign * Q[i, j])
    return m, mip


def unbounded_certificate(path):
    """Independent proof that a continuous relaxation is unbounded (or None if there is none): a
    feasible point of the rows and bounds and a recession direction d with Q d = 0 and c'd < 0, both
    from HiGHS's LP solver (|d| <= 1), checked again in floating point here."""
    import highspy

    m = mps.read_mps(path)
    n = m.num_vars
    c = np.array(m.obj, dtype=float) * m.sense
    Q = np.zeros((n, n))
    for (i, j), v in m.quad.items():
        Q[i, j] += v * m.sense
        if i != j:
            Q[j, i] += v * m.sense
    A = np.zeros((m.num_cons, n))
    for i in range(m.num_cons):
        for j, v in zip(m.row_cols[i], m.row_vals[i]):
            A[i, j] += v

    def solve_lp(cost, lb, ub, M, rlo, rhi):
        h = highspy.Highs()
        h.setOptionValue("output_flag", False)
        lp = highspy.HighsLp()
        lp.num_col_, lp.num_row_ = len(cost), M.shape[0]
        big = highspy.kHighsInf
        lp.col_cost_ = [float(v) for v in cost]
        lp.col_lower_ = [max(-big, v) for v in lb]
        lp.col_upper_ = [min(big, v) for v in ub]
        lp.row_lower_ = [max(-big, v) for v in rlo]
        lp.row_upper_ = [min(big, v) for v in rhi]
        lp.a_matrix_.format_ = highspy.MatrixFormat.kColwise
        starts, idx, val = [0], [], []
        for j in range(len(cost)):
            for i in np.nonzero(M[:, j])[0]:
                idx.append(int(i))
                val.append(float(M[i, j]))
            starts.append(len(idx))
        lp.a_matrix_.start_, lp.a_matrix_.index_, lp.a_matrix_.value_ = starts, idx, val
        h.passModel(lp)
        h.run()
        return h.modelStatusToString(h.getModelStatus()), np.array(h.getSolution().col_value)

    st, _ = solve_lp(np.zeros(n), m.lb, m.ub, A, m.row_lo, m.row_hi)
    if st != "Optimal":
        return None
    dl = [0.0 if np.isfinite(v) else -1.0 for v in m.lb]
    du = [0.0 if np.isfinite(v) else 1.0 for v in m.ub]
    rlo = [0.0 if np.isfinite(v) else -np.inf for v in m.row_lo] + [0.0] * n
    rhi = [0.0 if np.isfinite(v) else np.inf for v in m.row_hi] + [0.0] * n
    st, d = solve_lp(c, dl, du, np.vstack([A, Q]) if n else A, rlo, rhi)
    if st != "Optimal" or not (c @ d < -1e-6) or np.abs(Q @ d).max(initial=0.0) > 1e-9:
        return None
    return float(c @ d)


def run_popso(popso, path, sol, extra):
    out = subprocess.run([str(popso), str(path), "--threads", "1", "--time-limit", "30", "--gap", "1e-9",
                          "--mem-limit", "500", "--sol", str(sol)] + extra, capture_output=True, text=True,
                         timeout=120)
    st = re.search(r"^Status: (\S+)", out.stdout, re.M)
    ob = re.search(r"^Objective: (\S+)", out.stdout, re.M)
    if out.returncode != 0 or not st:
        return "ERROR " + (out.stderr or out.stdout)[-300:], None
    return st.group(1).lower(), (float(ob.group(1)) if ob and ob.group(1) != "-" else None)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--n", type=int, default=200)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--build", default=str(ROOT / "build"))
    ap.add_argument("--format", choices=["mps", "lp"], default="mps")
    ap.add_argument("--mip-frac", type=float, default=0.5, help="fraction of models with integer columns")
    ap.add_argument("--popso-args", default="", help="extra popso arguments")
    ap.add_argument("--out", help="directory for the models (default: temporary)")
    args = ap.parse_args()
    popso = Path(args.build) / "popso"
    out = Path(args.out) if args.out else Path(tempfile.mkdtemp(prefix="qp_fuzz_"))
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)
    extra = shlex.split(args.popso_args)
    bad = 0
    counts = {}
    for k in range(args.n):
        m, mip = random_model(rng, k, args.mip_frac)
        path = out / f"qpfuzz{k}.{args.format}"
        m.write(path)
        mps_path = out / f"qpfuzz{k}.mps"
        if args.format != "mps":
            m.write(mps_path)
        ref_name = "scip" if mip else "highs"
        ref = (refsolvers.solve_scip if mip else refsolvers.solve_highs)(str(mps_path), 30, 1)
        if not mip and ref["status"] not in ("optimal", "infeasible", "unbounded"):
            ref_name, ref = "scip", refsolvers.solve_scip(str(mps_path), 30, 1)
        rst, robj = ref["status"], ref["objective"]
        st, obj = run_popso(popso, path, out / "fuzz.sol", extra)
        if ref_name == "highs" and (st != rst or (rst == "optimal" and (obj is None or robj is None or
                                                                       abs(obj - robj) > 1e-6 * max(1.0, abs(robj))))):
            # HiGHS's QP solver has been seen to misreport maximization and unbounded QPs: SCIP decides.
            sref = refsolvers.solve_scip(str(mps_path), 30, 1)
            print(f"  qpfuzz{k}: HiGHS {rst} {robj} vs popso {st} {obj}; SCIP {sref['status']} {sref['objective']}")
            ref_name, rst, robj = "scip", sref["status"], sref["objective"]
        if st == "unbounded" and rst != "unbounded" and unbounded_certificate(mps_path) is not None:
            # The references are unreliable on unbounded QPs (HiGHS reports finite "optima", SCIP
            # runs into its time limit); an LP certificate decides.
            print(f"  qpfuzz{k}: {ref_name} {rst} {robj}, popso unbounded: certified by a recession direction")
            ref_name, rst, robj = "certificate", "unbounded", None
        if rst not in ("optimal", "infeasible", "unbounded", "inf_or_unbd"):
            rst = "undecided"
            if st.startswith("ERROR") or (st in ("optimal", "unbounded") and obj is None):
                bad += 1
                print(f"BAD qpfuzz{k}: popso {st} {obj}, {ref_name} undecided  ({path})")
                continue
        ok = rst == "undecided" or st == rst or \
            (rst == "inf_or_unbd" and st in ("infeasible", "unbounded", "inf_or_unbd")) or \
            (st == "inf_or_unbd" and rst in ("infeasible", "unbounded"))
        if ok and rst == "optimal":
            ok = obj is not None and robj is not None and abs(obj - robj) <= 1e-6 * max(1.0, abs(robj))
        if ok and obj is not None and st == "optimal":
            sol = solfile.read_sol(out / "fuzz.sol")
            v = mps.verify(mps.read_mps(mps_path), sol.values)
            ok = v.feasible and abs(v.objective - obj) <= 1e-7 * max(1.0, abs(obj))
            if not ok:
                print(f"  qpfuzz{k}: popso solution fails verification: {v.worst} (objective {v.objective})")
        key = f"{'miqp' if mip else 'qp'} {rst}"
        counts[key] = counts.get(key, 0) + 1
        if not ok:
            bad += 1
            print(f"BAD qpfuzz{k}: popso {st} {obj}, {ref_name} {rst} {robj}  ({path})")
        else:
            path.unlink()
            if mps_path.exists():
                mps_path.unlink()
    print(f"{args.n} models ({', '.join(f'{k} {v}' for k, v in sorted(counts.items()))}): "
          f"{bad} mismatch(es)" + (f"; kept in {out}" if bad else ""))
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
