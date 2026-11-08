"""Independent pure-Python MPS reader and solution verifier.

Used by reference.py (to compute instance statistics and to cross-check the
reference solvers' solutions) and by run.py (to verify our solver's output).
It deliberately does not depend on any solver library.

Supported: free and fixed MPS (whitespace-tokenized, so names must not contain
spaces), optional .gz, NAME, OBJSENSE (inline or as a section; MIN/MAX/
MINIMIZE/MAXIMIZE), OBJSENSE/OBJSENSE MAX, ROWS (N/L/G/E; first N row is the
objective, other N rows are dropped), COLUMNS with integer MARKER blocks, RHS
(RHS on the objective row = negated objective constant), RANGES, BOUNDS
(UP/LO/FX/FR/MI/PL/BV/LI/UI, SC/SI semi-continuous / semi-integer), SOS
(" S1 SOS name prio" or " S1 name" headers, "col weight" members),
INDICATORS (" IF row col value": the row only has to hold when the column
equals the value), quadratic objectives (QUADOBJ: one triangle, QMATRIX and
QSECTION [objective name]: the whole matrix; objective c'x + 1/2 x'Qx),
ENDATA.  Comment lines start with '*'.  Quadratic constraints (QCMATRIX,
QSECTION of another row) are rejected.

Semantics (common MPS conventions, as in HiGHS/SCIP):
  * Default bounds [0, +inf].  EXCEPTION (verified to match HiGHS 1.15 and
    SCIP 10): an integer column from a MARKER block that has
    no entry at all in the BOUNDS section is binary, [0, 1].  Any BOUNDS
    entry for it (even only LO or MI) restores the default ub = +inf.
  * UP with a negative value on a column whose lower bound was not set sets
    the lower bound to -inf (CPLEX convention).  NOTE: HiGHS 1.15 and
    SCIP keep lb = 0 instead (-> infeasible / read error); see
    generators/lp_sanity.py HAND_DERIVED_ONLY.
  * MI sets lb = -inf (ub unchanged); PL sets ub = +inf.
  * BV: integer, bounds [0, 1].  LI/UI: integer + lower/upper bound.
  * RANGES R on row with rhs b:  L: [b-|R|, b]; G: [b, b+|R|];
    E: R>0 -> [b, b+R], R<0 -> [b+R, b].
  * SC bound: semi-continuous (semi-integer on an integer column) with the
    value as upper bound: x = 0 or lb <= x <= ub.  SI: semi-integer.
  * SOS1: at most one member nonzero; SOS2: at most two, adjacent in the
    order of the weights.
"""

from __future__ import annotations

import gzip
import math
from dataclasses import dataclass, field

INF = math.inf


class MPSError(Exception):
    pass


@dataclass
class Model:
    name: str = ""
    sense: int = 1  # 1 = minimize, -1 = maximize
    obj_name: str = ""
    obj_const: float = 0.0
    col_names: list = field(default_factory=list)
    col_index: dict = field(default_factory=dict)
    obj: list = field(default_factory=list)
    lb: list = field(default_factory=list)
    ub: list = field(default_factory=list)
    is_int: list = field(default_factory=list)
    row_names: list = field(default_factory=list)
    row_index: dict = field(default_factory=dict)
    row_type: list = field(default_factory=list)  # 'L','G','E'
    row_lo: list = field(default_factory=list)
    row_hi: list = field(default_factory=list)
    # row-wise sparse matrix: list of (col indices list, values list)
    row_cols: list = field(default_factory=list)
    row_vals: list = field(default_factory=list)
    warnings: list = field(default_factory=list)
    semi: list = field(default_factory=list)  # per column: semi-continuous / semi-integer
    sos: list = field(default_factory=list)  # (name, type, [col], [weight])
    # (name, binary col, value, lo, hi, [col], [val])
    indicators: list = field(default_factory=list)
    # quadratic objective 1/2 x'Qx: {(i, j): Q_ij} with i >= j (lower triangle)
    quad: dict = field(default_factory=dict)
    # quadratic constraint terms: {row: {(i, j): c}} adding c * x_i * x_j (i >= j, no 1/2) to the row's activity
    qrows: dict = field(default_factory=dict)

    @property
    def num_vars(self):
        return len(self.col_names)

    @property
    def num_cons(self):
        return len(self.row_names)

    @property
    def num_nonzeros(self):
        return sum(len(c) for c in self.row_cols)

    @property
    def num_integers(self):
        return sum(1 for v in self.is_int if v)

    @property
    def num_binaries(self):
        return sum(
            1
            for j, v in enumerate(self.is_int)
            if v and self.lb[j] >= 0 and self.ub[j] <= 1
        )

    def stats(self):
        st = {
            "num_vars": self.num_vars,
            "num_cons": self.num_cons,
            "num_nonzeros": self.num_nonzeros,
            "num_integers": self.num_integers,
            "num_binaries": self.num_binaries,
            "sense": "max" if self.sense < 0 else "min",
        }
        if any(self.semi):
            st["num_semi"] = sum(1 for v in self.semi if v)
        if self.sos:
            st["num_sos"] = len(self.sos)
        if self.indicators:
            st["num_indicators"] = len(self.indicators)
        return st


def _open(path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="latin-1")
    return open(path, "rt", encoding="latin-1")


def _float(tok):
    try:
        return float(tok)
    except ValueError as e:
        raise MPSError(f"bad number {tok!r}") from e


def read_mps(path) -> Model:
    m = Model()
    section = None
    qrow = None  # the constraint row of a QCMATRIX section
    in_int = False
    rhs_vals = {}  # row name -> value
    range_vals = {}
    rows_entries = {}  # row name -> list[(col, val)]
    # bounds bookkeeping
    lb_set = []
    has_bound = []
    obj_rows_extra = set()
    col_obj = {}
    entries = []  # (row_name, col_idx, val)
    ind_specs = []  # (row name, col idx, value)

    def col_id(name, integer):
        idx = m.col_index.get(name)
        if idx is None:
            idx = len(m.col_names)
            m.col_index[name] = idx
            m.col_names.append(name)
            m.obj.append(0.0)
            m.lb.append(0.0)
            m.ub.append(INF)
            m.is_int.append(integer)
            m.semi.append(False)
            lb_set.append(False)
            has_bound.append(False)
        return idx

    with _open(path) as f:
        for lineno, raw in enumerate(f, 1):
            line = raw.rstrip("\n\r")
            # A comment has its '*' in column 1: irish-electricity's RHS lines start with "    *RHS*".
            if not line.strip() or line.startswith("*"):
                continue
            toks = line.split()
            if not line[0].isspace():
                # section header
                head = toks[0].upper()
                if head == "NAME":
                    m.name = toks[1] if len(toks) > 1 else ""
                    section = "NAME"
                elif head == "OBJSENSE":
                    section = "OBJSENSE"
                    if len(toks) > 1:
                        _set_sense(m, toks[1])
                elif head in ("ROWS", "COLUMNS", "RHS", "RANGES", "BOUNDS", "SOS", "INDICATORS"):
                    section = head
                elif head == "ENDATA":
                    break
                elif head in ("OBJNAME",):
                    section = "OBJNAME"
                    if len(toks) > 1:
                        m.obj_name = toks[1]
                elif head in ("QUADOBJ", "QMATRIX", "QSECTION"):
                    qrow = None
                    if head == "QSECTION" and len(toks) > 1 and toks[1] != m.obj_name:
                        if toks[1] not in m.row_index:
                            raise MPSError(f"line {lineno}: QSECTION for unknown row {toks[1]}")
                        qrow = m.row_index[toks[1]]
                    section = head if qrow is None else "QCMATRIX"
                elif head == "QCMATRIX":
                    if len(toks) != 2 or toks[1] not in m.row_index:
                        raise MPSError(f"line {lineno}: QCMATRIX needs a constraint row")
                    qrow = m.row_index[toks[1]]
                    section = "QCMATRIX"
                elif head in ("GENERAL", "PWLOBJ"):
                    raise MPSError(f"unsupported section {head}")
                else:
                    raise MPSError(f"line {lineno}: unknown section {head!r}")
                continue

            if section == "OBJSENSE":
                _set_sense(m, toks[0])
            elif section == "OBJNAME":
                m.obj_name = toks[0]
            elif section == "ROWS":
                t, name = toks[0].upper(), toks[1]
                if t == "N":
                    if not m.obj_name:
                        m.obj_name = name
                    elif name != m.obj_name:
                        obj_rows_extra.add(name)
                    continue
                if t not in ("L", "G", "E"):
                    raise MPSError(f"line {lineno}: bad row type {t}")
                if name in m.row_index:
                    raise MPSError(f"line {lineno}: duplicate row {name}")
                m.row_index[name] = len(m.row_names)
                m.row_names.append(name)
                m.row_type.append(t)
            elif section == "COLUMNS":
                if len(toks) >= 3 and toks[1].strip("'\"").upper() == "MARKER":
                    kind = toks[2].strip("'\"").upper()
                    if kind == "INTORG":
                        in_int = True
                    elif kind == "INTEND":
                        in_int = False
                    else:
                        raise MPSError(f"line {lineno}: bad marker {kind}")
                    continue
                cname = toks[0]
                j = col_id(cname, in_int)
                if in_int:
                    m.is_int[j] = True
                if len(toks) not in (1, 3, 5):
                    raise MPSError(f"line {lineno}: bad COLUMNS line")
                for k in range(1, len(toks) - 1, 2):
                    rname, val = toks[k], _float(toks[k + 1])
                    if rname == m.obj_name:
                        m.obj[j] += val
                    elif rname in obj_rows_extra:
                        pass
                    else:
                        if rname not in m.row_index:
                            raise MPSError(f"line {lineno}: unknown row {rname}")
                        entries.append((m.row_index[rname], j, val))
            elif section in ("RHS", "RANGES"):
                pairs = toks[1:] if len(toks) % 2 == 1 else toks
                target = rhs_vals if section == "RHS" else range_vals
                for k in range(0, len(pairs) - 1, 2):
                    rname, val = pairs[k], _float(pairs[k + 1])
                    if rname == m.obj_name:
                        if section == "RHS":
                            m.obj_const = -val
                        continue
                    if rname in obj_rows_extra:
                        continue
                    if rname not in m.row_index:
                        raise MPSError(f"line {lineno}: unknown row {rname}")
                    target[rname] = val
            elif section == "SOS":
                t0 = toks[0].upper()
                is_header = t0 in ("S1", "S2") and (len(toks) == 1 or toks[1].upper() == "SOS"
                                                     or not _is_number(toks[1]))
                if is_header:
                    k = 2 if len(toks) > 1 and toks[1].upper() == "SOS" else 1
                    sname = toks[k] if len(toks) > k else f"SOS{len(m.sos) + 1}"
                    m.sos.append((sname, int(t0[1]), [], []))
                    continue
                if not m.sos:
                    raise MPSError(f"line {lineno}: SOS member before a set header")
                if len(toks) == 1 and ":" in toks[0]:
                    cname, w = toks[0].split(":", 1)
                    w = _float(w)
                elif len(toks) == 2:
                    cname, w = toks[0], _float(toks[1])
                elif len(toks) == 1:
                    cname, w = toks[0], float(len(m.sos[-1][2]) + 1)
                else:
                    raise MPSError(f"line {lineno}: bad SOS line")
                if cname not in m.col_index:
                    raise MPSError(f"line {lineno}: unknown column {cname}")
                m.sos[-1][2].append(m.col_index[cname])
                m.sos[-1][3].append(w)
            elif section == "QCMATRIX":
                if len(toks) != 3:
                    raise MPSError(f"line {lineno}: bad QCMATRIX line")
                a, b = m.col_index.get(toks[0]), m.col_index.get(toks[1])
                if a is None or b is None:
                    raise MPSError(f"line {lineno}: unknown column in QCMATRIX")
                terms = m.qrows.setdefault(qrow, {})
                key = (max(a, b), min(a, b))
                terms[key] = terms.get(key, 0.0) + _float(toks[2])  # (a, b) and (b, a) both listed: they add up
            elif section in ("QUADOBJ", "QMATRIX", "QSECTION"):
                if len(toks) != 3:
                    raise MPSError(f"line {lineno}: bad {section} line")
                a, b = m.col_index.get(toks[0]), m.col_index.get(toks[1])
                if a is None or b is None:
                    raise MPSError(f"line {lineno}: unknown column in {section}")
                v = _float(toks[2])
                if section != "QUADOBJ" and a != b:
                    v *= 0.5  # both (a, b) and (b, a) are listed
                key = (max(a, b), min(a, b))
                m.quad[key] = m.quad.get(key, 0.0) + v
            elif section == "INDICATORS":
                if toks[0].upper() != "IF" or len(toks) not in (3, 4):
                    raise MPSError(f"line {lineno}: bad INDICATORS line")
                if toks[1] not in m.row_index or toks[2] not in m.col_index:
                    raise MPSError(f"line {lineno}: unknown row or column")
                val = int(_float(toks[3])) if len(toks) == 4 else 1
                ind_specs.append((toks[1], m.col_index[toks[2]], val))
            elif section == "BOUNDS":
                t = toks[0].upper()
                if t in ("SC", "SI"):
                    # semi-continuous / semi-integer; the value (if any) is the upper bound
                    if len(toks) == 2:
                        cname, val = toks[1], None
                    elif len(toks) == 3:
                        if toks[1] in m.col_index and _is_number(toks[2]):
                            cname, val = toks[1], _float(toks[2])
                        else:
                            cname, val = toks[2], None
                    else:
                        cname, val = toks[2], _float(toks[3])
                    if cname not in m.col_index:
                        raise MPSError(f"line {lineno}: unknown column {cname}")
                    j = m.col_index[cname]
                    has_bound[j] = True
                    m.semi[j] = True
                    if t == "SI":
                        m.is_int[j] = True
                    if val is not None:
                        m.ub[j] = INF if val >= 1e30 else val
                    continue
                novalue = t in ("FR", "MI", "PL", "BV")
                if novalue:
                    if len(toks) == 2:
                        cname, val = toks[1], None
                    elif len(toks) == 3:
                        # either "set col" or "col val"
                        if toks[2] in m.col_index or toks[1] not in m.col_index:
                            cname, val = toks[2], None
                        else:
                            cname, val = toks[1], None
                    else:
                        cname, val = toks[2], None
                else:
                    if len(toks) == 3:
                        cname, val = toks[1], _float(toks[2])
                    elif len(toks) >= 4:
                        cname, val = toks[2], _float(toks[3])
                    else:
                        raise MPSError(f"line {lineno}: bad BOUNDS line")
                if cname not in m.col_index:
                    raise MPSError(f"line {lineno}: unknown column {cname}")
                j = m.col_index[cname]
                has_bound[j] = True
                if t == "UP":
                    if val < 0 and m.lb[j] == 0 and not lb_set[j]:
                        m.lb[j] = -INF
                        m.warnings.append(f"negative UP bound on {cname}: lb set to -inf")
                    m.ub[j] = val
                elif t == "LO":
                    m.lb[j] = val
                    lb_set[j] = True
                elif t == "FX":
                    m.lb[j] = m.ub[j] = val
                    lb_set[j] = True
                elif t == "FR":
                    m.lb[j], m.ub[j] = -INF, INF
                    lb_set[j] = True
                elif t == "MI":
                    m.lb[j] = -INF
                    lb_set[j] = True
                elif t == "PL":
                    m.ub[j] = INF
                elif t == "BV":
                    m.is_int[j] = True
                    m.lb[j], m.ub[j] = 0.0, 1.0
                    lb_set[j] = True
                elif t == "LI":
                    m.is_int[j] = True
                    m.lb[j] = val
                    lb_set[j] = True
                elif t == "UI":
                    m.is_int[j] = True
                    if val < 0 and m.lb[j] == 0 and not lb_set[j]:
                        m.lb[j] = -INF
                    m.ub[j] = val
                else:
                    raise MPSError(f"line {lineno}: unsupported bound type {t}")
            else:
                raise MPSError(f"line {lineno}: data outside of a section")

    for j in range(len(m.col_names)):
        if m.is_int[j] and not has_bound[j]:
            m.ub[j] = 1.0  # MARKER integer without any bound entry: binary
    nrows = len(m.row_names)
    m.row_cols = [[] for _ in range(nrows)]
    m.row_vals = [[] for _ in range(nrows)]
    for i, j, v in entries:
        m.row_cols[i].append(j)
        m.row_vals[i].append(v)
    for i, name in enumerate(m.row_names):
        b = rhs_vals.get(name, 0.0)
        t = m.row_type[i]
        lo, hi = {"L": (-INF, b), "G": (b, INF), "E": (b, b)}[t]
        if name in range_vals:
            r = range_vals[name]
            if t == "L":
                lo = b - abs(r)
            elif t == "G":
                hi = b + abs(r)
            else:
                if r >= 0:
                    hi = b + r
                else:
                    lo = b + r
        m.row_lo.append(lo)
        m.row_hi.append(hi)
    if ind_specs:
        # indicator rows are not ordinary rows
        drop = set()
        for rname, b, val in ind_specs:
            i = m.row_index[rname]
            drop.add(i)
            m.indicators.append((rname, b, val, m.row_lo[i], m.row_hi[i], m.row_cols[i], m.row_vals[i]))
        keep = [i for i in range(nrows) if i not in drop]
        m.row_names = [m.row_names[i] for i in keep]
        m.row_type = [m.row_type[i] for i in keep]
        m.row_lo = [m.row_lo[i] for i in keep]
        m.row_hi = [m.row_hi[i] for i in keep]
        m.row_cols = [m.row_cols[i] for i in keep]
        m.row_vals = [m.row_vals[i] for i in keep]
        m.row_index = {n: k for k, n in enumerate(m.row_names)}
    return m


def _is_number(tok):
    try:
        float(tok)
        return True
    except ValueError:
        return False


def _set_sense(m, tok):
    t = tok.upper()
    if t in ("MAX", "MAXIMIZE", "MAXIMISE"):
        m.sense = -1
    elif t in ("MIN", "MINIMIZE", "MINIMISE"):
        m.sense = 1
    else:
        raise MPSError(f"bad OBJSENSE {tok!r}")


def objective(m: Model, x) -> float:
    lin = math.fsum(c * x[j] for j, c in enumerate(m.obj) if c)
    quad = math.fsum((0.5 if i == j else 1.0) * v * x[i] * x[j] for (i, j), v in m.quad.items())
    return m.obj_const + lin + quad


@dataclass
class Verification:
    feasible: bool
    objective: float
    max_bound_viol: float
    max_int_viol: float
    max_row_viol: float
    worst: str  # description of the worst violation
    missing: int  # number of variables missing from the solution
    max_logic_viol: float = 0.0  # SOS sets and active indicator rows

    def to_dict(self):
        return dict(self.__dict__)


def _scaled_viol(viol, bound):
    return viol / max(1.0, abs(bound))


def _cone_violation(m, i, x):
    """Row i as a second-order cone over its columns' bounds (as popso's solve_oa detect_cone): in its <= direction
    sum_k c_k (x_k + d_k)^2 + gamma <= c_t (t - s)^2 (t - s of one sign over t's bounds) or <= c u v (u, v >= 0),
    c > 0, gamma >= 0, no other linear terms. Returns its norm-form violation (sqrt(gamma + sum ...) - head) relative
    to the row's constants (max(1, |side|, shifts)) as a linear row's to its side, or None if the row is no such cone.
    The squared form's tolerance would let a column reach sqrt(tol) at the apex (x^2 <= t^2 + 1e-9 at t = 0:
    |x| = 3e-5); relative to the head, a rotated cone written as (u - v)^2 + 4 x^2 <= (u + v)^2 lets x reach 1.5 at
    u = 0, v = 6575 (CBLIB clay0304h)."""
    lo, hi = m.row_lo[i], m.row_hi[i]
    if math.isfinite(lo) and math.isfinite(hi):
        return None
    sgn = 1.0 if math.isfinite(hi) else -1.0
    side = sgn * (hi if sgn > 0 else lo)
    pos, t_col, uv = {}, None, None
    for (a, b), c in m.qrows[i].items():
        c *= sgn
        if c == 0.0:
            continue
        if a == b:
            if c > 0:
                pos[a] = pos.get(a, 0.0) + c
            elif t_col is None:
                t_col, ct = a, -c
            else:
                return None
        elif c < 0 and uv is None:
            uv, cuv = (a, b), -c
        else:
            return None
    if (t_col is None) == (uv is None):
        return None
    lin = {}
    for j, v in zip(m.row_cols[i], m.row_vals[i]):
        lin[j] = lin.get(j, 0.0) + sgn * v
    H, r2, scale = side, 0.0, max(1.0, abs(side))
    for j, c in pos.items():
        if uv is not None and j in uv:
            return None
        bj = lin.pop(j, 0.0)
        H += bj * bj / (4.0 * c)
        r2 += c * (x[j] + bj / (2.0 * c)) ** 2
        scale = max(scale, abs(bj / (2.0 * c)))
    if t_col is not None:
        bt = lin.pop(t_col, 0.0)
        sh = bt / (2.0 * ct)
        H -= bt * bt / (4.0 * ct)
        scale = max(scale, abs(sh))
        if m.lb[t_col] >= sh:
            head = math.sqrt(ct) * (x[t_col] - sh)
        elif m.ub[t_col] <= sh:
            head = math.sqrt(ct) * (sh - x[t_col])
        else:
            return None
    else:
        u, v = uv
        if not (m.lb[u] >= 0 and m.lb[v] >= 0):
            return None
        r2 += 0.25 * cuv * (x[u] - x[v]) ** 2
        head = 0.5 * math.sqrt(cuv) * (x[u] + x[v])
    if any(c != 0.0 for c in lin.values()) or H > 1e-12 * max(1.0, abs(side)):
        return None
    return (math.sqrt(max(0.0, r2 - H)) - head) / scale


def verify(m: Model, values: dict, feastol=1e-6, inttol=1e-6) -> Verification:
    """Check a solution given as {var_name: value}.

    Violations are measured relative to max(1, |bound|) (i.e. absolute for
    small bounds, relative for large ones).  Unknown names are ignored;
    missing variables are taken as 0 and counted in `missing`.
    """
    x = [0.0] * m.num_vars
    missing = 0
    for j, name in enumerate(m.col_names):
        v = values.get(name)
        if v is None:
            missing += 1
        else:
            x[j] = v
    worst = ""
    worst_v = 0.0
    mb = mi = mr = 0.0
    for j in range(m.num_vars):
        xj = x[j]
        if not math.isfinite(xj):
            return Verification(False, math.nan, INF, INF, INF,
                                f"non-finite value for {m.col_names[j]}", missing)
        if m.semi[j] and abs(xj) <= feastol:
            pass  # semi-continuous at 0
        elif xj < m.lb[j]:
            v = _scaled_viol(m.lb[j] - xj, m.lb[j])
            if v > mb:
                mb = v
            if v > worst_v:
                worst_v, worst = v, f"var {m.col_names[j]}={xj} < lb {m.lb[j]}"
        if xj > m.ub[j]:
            v = _scaled_viol(xj - m.ub[j], m.ub[j])
            if v > mb:
                mb = v
            if v > worst_v:
                worst_v, worst = v, f"var {m.col_names[j]}={xj} > ub {m.ub[j]}"
        if m.is_int[j]:
            v = abs(xj - round(xj))
            if v > mi:
                mi = v
            if v > inttol and v > worst_v:
                worst_v, worst = v, f"var {m.col_names[j]}={xj} not integral"
    for i in range(m.num_cons):
        if i in m.qrows:
            cv = _cone_violation(m, i, x)
            if cv is not None:
                if cv > mr:
                    mr = cv
                if cv > worst_v:
                    worst_v, worst = cv, f"cone row {m.row_names[i]} violated by {cv:.3g} (norm form, relative)"
                continue
        act = math.fsum(v * x[j] for j, v in zip(m.row_cols[i], m.row_vals[i]))
        lo, hi = m.row_lo[i], m.row_hi[i]
        # Floating point: terms of 5e9 summing to 0 carry ~1e-6 of rounding (QPLIB_9002: a solution's row at 1.0003e-6),
        # so a row's scale includes 1e-9 of its terms' magnitude (no effect below terms of ~1e9).
        lmag = 1e-9 * math.fsum(abs(v * x[j]) for j, v in zip(m.row_cols[i], m.row_vals[i]))
        scale_lo = max(abs(lo), lmag) if math.isfinite(lo) else lo
        scale_hi = max(abs(hi), lmag) if math.isfinite(hi) else hi
        if i in m.qrows:
            # A quadratic row: relative to its terms' magnitude too (they may cancel).
            q = [c * x[a] * x[b] for (a, b), c in m.qrows[i].items()]
            act += math.fsum(q)
            mag = math.fsum(abs(v * x[j]) for j, v in zip(m.row_cols[i], m.row_vals[i])) + math.fsum(abs(t) for t in q)
            scale_lo = max(abs(lo) if math.isfinite(lo) else 0.0, mag)
            scale_hi = max(abs(hi) if math.isfinite(hi) else 0.0, mag)
        if act < lo:
            v = _scaled_viol(lo - act, scale_lo)
            if v > mr:
                mr = v
            if v > worst_v:
                worst_v, worst = v, f"row {m.row_names[i]} activity {act} < {lo}"
        if act > hi:
            v = _scaled_viol(act - hi, scale_hi)
            if v > mr:
                mr = v
            if v > worst_v:
                worst_v, worst = v, f"row {m.row_names[i]} activity {act} > {hi}"
    ml = 0.0
    for sname, t, cols, weights in m.sos:
        vals = [abs(x[j]) for _, j in sorted(zip(weights, cols))]
        if len(vals) <= t:
            continue
        if t == 1:
            v = sorted(vals)[-2]
        else:  # largest member outside the best pair of neighbours
            v = min(max(vals[:k] + vals[k + 2:], default=0.0) for k in range(len(vals) - 1))
        if v > ml:
            ml = v
        if v > worst_v:
            worst_v, worst = v, f"SOS{t} {sname} violated by {v}"
    for iname, b, val, lo, hi, cols, vals in m.indicators:
        if abs(x[b] - val) >= 0.5:
            continue
        act = math.fsum(a * x[j] for j, a in zip(cols, vals))
        v = max(_scaled_viol(lo - act, lo) if act < lo else 0.0, _scaled_viol(act - hi, hi) if act > hi else 0.0)
        if v > ml:
            ml = v
        if v > worst_v:
            worst_v, worst = v, f"indicator {iname} activity {act} outside [{lo}, {hi}]"
    feasible = mb <= feastol and mr <= feastol and mi <= inttol and ml <= feastol
    return Verification(feasible, objective(m, x), mb, mi, mr, worst, missing, ml)
