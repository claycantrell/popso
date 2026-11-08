#!/usr/bin/env python3
"""Download mixed-integer second-order cone instances of CBLIB (cblib.zib.de) and convert them to MPS with
quadratic constraints (QCMATRIX) for the benchmark harness.

    .venv/bin/python benchmarks/tools/cblib.py [--out benchmarks/cblib] [--check] [NAMES or name prefixes ...]

Without names a default selection of MISOCP families is fetched (a few instances of each portfolio family, the
others whole). CBF (the conic benchmark format, versions 1-3) puts the variables x in cones (F free, L+ / L- / L=
orthants, Q quadratic, QR rotated quadratic) and the affine rows A x + b in cones. Instances with semidefinite,
exponential or power cones, or without integer variables, are skipped. The MPS model:
- a variable cone L+ / L- / L= gives bounds; a quadratic cone (x1, ..., xk) the row x2^2 + ... + xk^2 - x1^2 <= 0
  with x1 >= 0, a rotated one (2 x1 x2 >= x3^2 + ... + xk^2) the row x3^2 + ... + xk^2 - 2 x1 x2 <= 0 with
  x1, x2 >= 0 (popso and SCIP both recognize these as second-order cones);
- a row cone L+ / L- / L= gives a linear row (one with a single column becomes that column's bound); a quadratic
  cone on rows gives one column per row (the row's column itself when the row is a single column with coefficient
  1 and no constant, else a new free column y = a x + b defined by an equality row) and the cone's row over these.
Writes <out>/<name>.mps.gz, <out>/meta.json and benchmarks/sets/cblib_misocp.txt. With --check every converted
instance is also solved by SCIP from the CBF file and from the MPS file (--check-time each) and the two results
compared (a conversion check independent of popso).
"""

from __future__ import annotations

import argparse
import gzip
import json
import math
import re
import sys
import tempfile
import urllib.request
from pathlib import Path

BENCH = Path(__file__).resolve().parent.parent
REPO = BENCH.parent
sys.path.insert(0, str(BENCH))

from lib.modeler import Model, lin, quicksum  # noqa: E402

SITE = "https://cblib.zib.de/download/all/"
INF = math.inf
CONES = {"F", "L+", "L-", "L=", "Q", "QR"}

# Default selection: name prefixes (all instances matching), and families of many instances (the first few).
DEFAULT_PREFIXES = ["clay", "flay", "slay", "estein", "pp-n", "sssd-", "uflquad-", "tls", "turbine", "nvs", "enpro",
                    "batch", "fac", "du-opt", "m3", "m6", "m7", "o7", "o8", "o9", "no7", "port_", "achtziger_stolpe",
                    "b1bigflowc", "netmod", "isil", "jha", "gams", "ex1223", "fiac", "beck", "demb", "gbd", "car",
                    "fang", "cx02", "qssp", "sched", "squfl", "ravem", "sambal", "varun", "syn", "rsyn"]
DEFAULT_FAMILIES = {"classical_50_": 5, "robust_50_": 5, "shortfall_50_": 5, "classical_200_": 2, "robust_100_": 2,
                    "shortfall_100_": 2}


class Unsupported(Exception):
    pass


def parse_cbf(text):
    """The CBF sections popso can take, as a dict; raises Unsupported for the others."""
    lines = [ln.strip() for ln in text.splitlines()]
    lines = [ln for ln in lines if ln and not ln.startswith("#")]
    d = dict(sense="MIN", n=0, m=0, var=[], con=[], ints=[], obja={}, objb=0.0, a=[], b={})
    k = 0

    def nxt():
        nonlocal k
        k += 1
        return lines[k - 1]

    while k < len(lines):
        kw = nxt()
        if kw == "VER":
            nxt()
        elif kw == "OBJSENSE":
            d["sense"] = nxt().split()[0]
        elif kw in ("VAR", "CON"):
            size, chunks = map(int, nxt().split())
            d["n" if kw == "VAR" else "m"] = size
            for _ in range(chunks):
                cone, sz = nxt().split()
                if cone not in CONES:
                    raise Unsupported(f"cone {cone}")
                d["var" if kw == "VAR" else "con"].append((cone, int(sz)))
        elif kw == "INT":
            d["ints"] = [int(nxt()) for _ in range(int(nxt()))]
        elif kw == "OBJACOORD":
            for _ in range(int(nxt())):
                j, a = nxt().split()
                d["obja"][int(j)] = d["obja"].get(int(j), 0.0) + float(a)
        elif kw == "OBJBCOORD":
            d["objb"] = float(nxt())
        elif kw == "ACOORD":
            for _ in range(int(nxt())):
                i, j, a = nxt().split()
                d["a"].append((int(i), int(j), float(a)))
        elif kw == "BCOORD":
            for _ in range(int(nxt())):
                i, b = nxt().split()
                d["b"][int(i)] = d["b"].get(int(i), 0.0) + float(b)
        else:
            raise Unsupported(f"section {kw}")
    if not d["ints"]:
        raise Unsupported("no integer variables")
    return d


def build(d, name):
    md = Model(name)
    md.sense = "max" if d["sense"] == "MAX" else "min"
    md.obj_const = d["objb"]
    n, m = d["n"], d["m"]
    lb, ub = [-INF] * n, [INF] * n
    var_cones = []
    j = 0
    for cone, sz in d["var"]:
        idx = list(range(j, j + sz))
        j += sz
        if cone == "L+":
            for i in idx:
                lb[i] = 0.0
        elif cone == "L-":
            for i in idx:
                ub[i] = 0.0
        elif cone == "L=":
            for i in idx:
                lb[i] = ub[i] = 0.0
        elif cone in ("Q", "QR"):
            for i in idx[: 1 if cone == "Q" else 2]:
                lb[i] = max(lb[i], 0.0)
            var_cones.append((cone, idx))
    rows = [dict() for _ in range(m)]
    for i, jj, a in d["a"]:
        rows[i][jj] = rows[i].get(jj, 0.0) + a
    b = d["b"]
    lin_rows, row_cones = [], []
    i0 = 0
    for cone, sz in d["con"]:
        idx = list(range(i0, i0 + sz))
        i0 += sz
        if cone == "F":
            continue
        if cone in ("Q", "QR"):
            row_cones.append((cone, idx))
            continue
        for i in idx:
            terms = {jj: a for jj, a in rows[i].items() if a != 0.0}
            c = b.get(i, 0.0)
            if len(terms) == 1:
                ((jj, a),) = terms.items()
                v = -c / a  # a x + c (cone) 0  <=>  x (cone, flipped for a < 0) v
                lower = (cone == "L+") == (a > 0)
                if cone == "L=":
                    lb[jj], ub[jj] = max(lb[jj], v), min(ub[jj], v)
                elif lower:
                    lb[jj] = max(lb[jj], v)
                else:
                    ub[jj] = min(ub[jj], v)
            elif not terms:
                ok = {"L+": c >= -1e-9, "L-": c <= 1e-9, "L=": abs(c) <= 1e-9}[cone]
                if not ok:
                    raise Unsupported(f"row {i} is an infeasible constant")
            else:
                lin_rows.append((i, cone, terms, c))
    ints = set(d["ints"])
    xs = []
    for jj in range(n):
        lo, hi = lb[jj], ub[jj]
        if jj in ints:
            lo, hi = math.ceil(lo - 1e-9) if lo > -INF else lo, math.floor(hi + 1e-9) if hi < INF else hi
        xs.append(md.add_var(lb=lo, ub=hi, obj=d["obja"].get(jj, 0.0), vtype="I" if jj in ints else "C", name=f"x{jj}"))
    for i, cone, terms, c in lin_rows:
        e = lin(list(terms.values()), [xs[jj] for jj in terms])
        md.add_constr(e >= -c if cone == "L+" else e <= -c if cone == "L-" else e == -c, f"r{i}")
    ncone = 0

    def cone_row(cone, members):
        nonlocal ncone
        heads, tails = (members[:1], members[1:]) if cone == "Q" else (members[:2], members[2:])
        if not tails:
            return  # only the heads' sign
        con = md.add_constr(quicksum([]) <= 0.0, f"k{ncone}")
        ncone += 1
        for t in tails:
            md.add_qterm(con, t, t, 1.0)
        if cone == "Q":
            md.add_qterm(con, heads[0], heads[0], -1.0)
        else:
            md.add_qterm(con, heads[0], heads[1], -2.0)

    for cone, idx in var_cones:
        cone_row(cone, [xs[jj] for jj in idx])
    for cone, idx in row_cones:
        members = []
        nheads = 1 if cone == "Q" else 2
        for pos, i in enumerate(idx):
            terms = {jj: a for jj, a in rows[i].items() if a != 0.0}
            c = b.get(i, 0.0)
            if c == 0.0 and len(terms) == 1 and next(iter(terms.values())) == 1.0:
                v = xs[next(iter(terms))]
                if pos < nheads:
                    v.lb = max(v.lb, 0.0)
            else:
                v = md.add_var(lb=0.0 if pos < nheads else -INF, ub=INF, name=f"y{i}")
                md.add_constr(v - lin(list(terms.values()), [xs[jj] for jj in terms]) == c, f"d{i}")
            members.append(v)
        cone_row(cone, members)
    return md


def scip_solve(path, time_limit):
    import pyscipopt

    s = pyscipopt.Model()
    s.hideOutput()
    s.readProblem(str(path))
    s.setParam("limits/time", time_limit)
    s.optimize()
    return s.getStatus(), (s.getObjVal() if s.getNSols() else None), s.getDualbound()


def listing():
    t = urllib.request.urlopen(SITE, timeout=120).read().decode("utf-8", "replace")
    return re.findall(r'<a href="([^"]+)\.cbf\.gz">', t)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default=str(BENCH / "cblib"))
    ap.add_argument("--set", default=str(BENCH / "sets" / "cblib_misocp.txt"), help="set file to write ('' = none)")
    ap.add_argument("--check", action="store_true", help="compare SCIP on the CBF and on the MPS file")
    ap.add_argument("--check-time", type=float, default=60.0)
    ap.add_argument("--max-cols", type=int, default=100000, help="skip instances with more variables")
    ap.add_argument("names", nargs="*", help="instance names or name prefixes (default: the MISOCP selection)")
    args = ap.parse_args(argv)

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    meta_path = out / "meta.json"
    meta = json.loads(meta_path.read_text()) if meta_path.exists() else {}
    names = listing()
    if args.names:
        todo = [nm for nm in names if any(nm == p or nm.startswith(p) for p in args.names)]
    else:
        todo = [nm for nm in names if any(nm.startswith(p) for p in DEFAULT_PREFIXES)]
        for fam, count in DEFAULT_FAMILIES.items():
            todo += sorted((nm for nm in names if nm.startswith(fam)), key=lambda s: int(s[len(fam):]))[:count]
    print(f"{len(todo)} instances selected", flush=True)
    with tempfile.TemporaryDirectory() as tmp:
        for nm in todo:
            rec = dict(name=nm)
            try:
                raw = urllib.request.urlopen(f"{SITE}{nm}.cbf.gz", timeout=300).read()
                d = parse_cbf(gzip.decompress(raw).decode("utf-8", "replace"))
                if d["n"] > args.max_cols:
                    raise Unsupported(f"{d['n']} variables (--max-cols {args.max_cols})")
                md = build(d, nm)
                text = md.to_mps()
                (out / f"{nm}.mps.gz").write_bytes(gzip.compress(text.encode()))
                rec.update(ok=True, vars=len(md.vars), rows=len(md.constrs),
                           ints=sum(v.vtype == "I" for v in md.vars),
                           cones=sum(1 for c in md.constrs if c.quad))
                if args.check:  # needs a SCIP with its CBF reader
                    try:
                        cbf = Path(tmp) / f"{nm}.cbf"
                        cbf.write_bytes(gzip.decompress(raw))
                        mps = Path(tmp) / f"{nm}.mps"
                        mps.write_text(text)
                        a, b_ = scip_solve(cbf, args.check_time), scip_solve(mps, args.check_time)
                        rec.update(scip_cbf=a, scip_mps=b_)
                        if a[0] == "optimal" and b_[0] == "optimal":
                            rec["check"] = abs(a[1] - b_[1]) <= 1e-5 * max(1.0, abs(a[1]))
                        else:
                            rec["check"] = None
                    except Exception as e:
                        rec.update(check=None, check_error=str(e)[:120])
            except Unsupported as e:
                rec.update(ok=False, skip=str(e))
            except Exception as e:  # download or parse failure
                rec.update(ok=False, error=f"{type(e).__name__}: {e}"[:200])
            meta[nm] = rec
            meta_path.write_text(json.dumps(meta, indent=1, sort_keys=True) + "\n")
            print(nm, {k: v for k, v in rec.items() if k != "name"}, flush=True)
    if args.set:
        lines = ["# CBLIB mixed-integer SOCP instances (benchmarks/tools/cblib.py): name, columns, rows, cones"]
        rel = out.resolve()
        for nm, r in sorted(meta.items()):
            if r.get("ok") and r.get("check") is not False:
                p = rel / f"{nm}.mps.gz"
                lines.append(f"{p.relative_to(REPO) if p.is_relative_to(REPO) else p}  # {r['vars']} {r['rows']} {r['cones']}")
        Path(args.set).write_text("\n".join(lines) + "\n")
        print(f"wrote {args.set}: {len(lines) - 1} instances")
    return 0


if __name__ == "__main__":
    sys.exit(main())
