#!/usr/bin/env python3
"""Download the Maros-Meszaros convex QP test set (the .mat files that OSQP's benchmarks and
qpsolvers/maros_meszaros_qpbenchmark distribute) and convert them to free MPS with a QUADOBJ
section, for testing popso's QP solver.

    python3 -m venv /tmp/mmvenv && /tmp/mmvenv/bin/pip install scipy numpy
    /tmp/mmvenv/bin/python benchmarks/tools/maros_meszaros.py --out /tmp/maros [--max-kb 60] [NAMES ...]

The .mat files hold  min 1/2 x'Px + q'x + r  s.t.  l <= Ax <= u  (|l|, |u| >= 1e19: infinite);
rows with one entry become column bounds, the others ranged / one-sided / equality rows; P's
lower triangle is written as QUADOBJ. Needs scipy (not used anywhere else in the repository,
hence the separate virtual environment). tests/data/maros_meszaros/ holds a few converted
problems used by the unit tests.
"""

from __future__ import annotations

import argparse
import json
import sys
import urllib.request
from pathlib import Path

REPO = "qpsolvers/maros_meszaros_qpbenchmark"
RAW = f"https://raw.githubusercontent.com/{REPO}/main/data/"
API = f"https://api.github.com/repos/{REPO}/contents/data"
BIG = 1e19


def listing():
    with urllib.request.urlopen(API) as f:
        return [(e["name"], e["size"]) for e in json.load(f) if e["name"].endswith(".mat")]


def fmt(v):
    return repr(float(v))


def convert(src: Path, dst: Path):
    import numpy as np
    import scipy.io
    import scipy.sparse as sp

    d = scipy.io.loadmat(src)
    P = sp.csc_matrix(d["P"])
    q = np.asarray(d["q"], dtype=float).ravel()
    r = float(np.asarray(d["r"]).ravel()[0]) if "r" in d else 0.0
    A = sp.csr_matrix(d["A"])
    lo = np.asarray(d["l"], dtype=float).ravel()
    up = np.asarray(d["u"], dtype=float).ravel()
    n = P.shape[0]
    clb, cub = np.full(n, -np.inf), np.full(n, np.inf)
    rows = []
    for i in range(A.shape[0]):
        s, e = A.indptr[i], A.indptr[i + 1]
        idx, val = A.indices[s:e], A.data[s:e]
        li = -np.inf if lo[i] <= -BIG else lo[i]
        ui = np.inf if up[i] >= BIG else up[i]
        if len(idx) == 1 and val[0] != 0:
            j, a = idx[0], val[0]
            bl, bu = (li / a, ui / a) if a > 0 else (ui / a, li / a)
            clb[j], cub[j] = max(clb[j], bl), min(cub[j], bu)
        elif len(idx) > 0:
            rows.append((idx, val, li, ui))
    cols = [[] for _ in range(n)]
    for k, (idx, val, _, _) in enumerate(rows):
        for j, a in zip(idx, val):
            cols[j].append((k, a))
    out = [f"NAME {src.stem}", "ROWS", " N  OBJ"]
    for k, (_, _, li, ui) in enumerate(rows):
        t = "E" if li == ui else ("N" if li == -np.inf and ui == np.inf else "L" if li == -np.inf else "G")
        out.append(f" {t}  R{k}")
    out.append("COLUMNS")
    for j in range(n):
        out.append(f"    C{j}  OBJ  {fmt(q[j])}")
        out.extend(f"    C{j}  R{k}  {fmt(a)}" for k, a in cols[j])
    out.append("RHS")
    if r != 0.0:
        out.append(f"    RHS  OBJ  {fmt(-r)}")
    ranges = []
    for k, (_, _, li, ui) in enumerate(rows):
        if li == ui or (li == -np.inf and ui != np.inf):
            out.append(f"    RHS  R{k}  {fmt(ui)}")
        elif li != -np.inf:
            out.append(f"    RHS  R{k}  {fmt(li)}")
            if ui != np.inf:
                ranges.append(f"    RNG  R{k}  {fmt(ui - li)}")
    if ranges:
        out += ["RANGES"] + ranges
    out.append("BOUNDS")
    for j in range(n):
        lb, ub = clb[j], cub[j]
        if lb == ub:
            out.append(f" FX BND  C{j}  {fmt(lb)}")
        elif lb == -np.inf and ub == np.inf:
            out.append(f" FR BND  C{j}")
        else:
            out.append(f" MI BND  C{j}" if lb == -np.inf else f" LO BND  C{j}  {fmt(lb)}")
            if ub != np.inf:
                out.append(f" UP BND  C{j}  {fmt(ub)}")
    out.append("QUADOBJ")
    L = sp.tril(P, format="csc")
    for j in range(n):
        for p in range(L.indptr[j], L.indptr[j + 1]):
            if L.data[p] != 0:
                out.append(f"    C{j}  C{L.indices[p]}  {fmt(L.data[p])}")
    out.append("ENDATA")
    dst.write_text("\n".join(out) + "\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("names", nargs="*", help="problem names (default: all up to --max-kb)")
    ap.add_argument("--out", required=True, help="output directory (MPS files; the .mat files go to <out>/mat)")
    ap.add_argument("--max-kb", type=float, default=60.0, help="skip .mat files larger than this (default 60)")
    args = ap.parse_args()
    out = Path(args.out)
    (out / "mat").mkdir(parents=True, exist_ok=True)
    names = args.names or [n[:-4] for n, size in listing() if size <= args.max_kb * 1024]
    for name in names:
        mat = out / "mat" / f"{name}.mat"
        if not mat.exists():
            urllib.request.urlretrieve(RAW + f"{name}.mat", mat)
        try:
            convert(mat, out / f"{name}.mps")
        except Exception as e:  # noqa: BLE001 (report and go on)
            print(f"{name}: conversion failed: {e}", file=sys.stderr)
    print(f"{len(names)} problems in {out}")


if __name__ == "__main__":
    main()
