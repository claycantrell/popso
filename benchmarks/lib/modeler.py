"""Tiny pure-Python linear modelling layer with an MPS writer.

Just enough to build the random benchmark families (generators/families.py)
without any solver library:

    m = Model("knap")
    x = [m.add_var(vtype="B", obj=p[j], name=f"x{j}") for j in range(n)]
    m.add_constr(quicksum(w[j] * x[j] for j in range(n)) <= C, "cap")
    m.sense = "max"
    m.write("knap.mps")

Variables have bounds (lb, ub; +-inf allowed), an objective coefficient and
a type: "C" continuous, "I" integer, "B" binary (integer with bounds [0, 1]).
Linear expressions support +, -, scalar * and quicksum(); comparing an
expression with <=, >= or == gives a constraint for Model.add_constr.  All
terms are moved to the left-hand side and constants to the right-hand side;
duplicate variables in one constraint are merged.

The MPS writer emits fixed-column-style free MPS: integer columns inside
MARKER INTORG/INTEND blocks, every integer column with explicit bounds
(BV for binaries, LI for an integer with an infinite upper bound), so the
file reads the same under all common MPS conventions.
"""

from __future__ import annotations

import math

INF = math.inf


class Var:
    __slots__ = ("index", "name", "lb", "ub", "obj", "vtype")
    __array_ufunc__ = None  # make numpy scalars defer to our operators

    def __init__(self, index, name, lb, ub, obj, vtype):
        self.index = index
        self.name = name
        self.lb = lb
        self.ub = ub
        self.obj = obj
        self.vtype = vtype

    def _expr(self):
        return LinExpr({self: 1.0})

    __hash__ = object.__hash__

    def __add__(self, o):
        return self._expr() + o

    __radd__ = __add__

    def __sub__(self, o):
        return self._expr() - o

    def __rsub__(self, o):
        return o + (-1.0) * self._expr()

    def __mul__(self, c):
        return self._expr() * c

    __rmul__ = __mul__

    def __neg__(self):
        return self._expr() * -1.0

    def __le__(self, o):
        return self._expr() <= o

    def __ge__(self, o):
        return self._expr() >= o

    def __eq__(self, o):
        return self._expr() == o

    def __repr__(self):
        return f"Var({self.name})"


class LinExpr:
    __slots__ = ("terms", "const")
    __array_ufunc__ = None

    def __init__(self, terms=None, const=0.0):
        self.terms = dict(terms) if terms else {}  # Var -> coefficient
        self.const = float(const)

    @staticmethod
    def _of(o):
        if isinstance(o, LinExpr):
            return o
        if isinstance(o, Var):
            return o._expr()
        return LinExpr(const=float(o))

    def __add__(self, o):
        o = LinExpr._of(o)
        t = dict(self.terms)
        for v, c in o.terms.items():
            t[v] = t.get(v, 0.0) + c
        return LinExpr(t, self.const + o.const)

    __radd__ = __add__

    def __sub__(self, o):
        return self + LinExpr._of(o) * -1.0

    def __rsub__(self, o):
        return LinExpr._of(o) + self * -1.0

    def __mul__(self, c):
        if isinstance(c, (Var, LinExpr)):
            raise TypeError("only linear expressions are supported")
        c = float(c)
        return LinExpr({v: c * a for v, a in self.terms.items()}, c * self.const)

    __rmul__ = __mul__

    def __neg__(self):
        return self * -1.0

    def __le__(self, o):
        return Constraint(self - o, "L")

    def __ge__(self, o):
        return Constraint(self - o, "G")

    def __eq__(self, o):
        return Constraint(self - o, "E")

    __hash__ = None


def quicksum(items):
    t, const = {}, 0.0
    for it in items:
        it = LinExpr._of(it)
        for v, c in it.terms.items():
            t[v] = t.get(v, 0.0) + c
        const += it.const
    return LinExpr(t, const)


def lin(coefs, vars_):
    """sum_j coefs[j] * vars_[j] (coefficients converted to float)."""
    return LinExpr(_merge(zip(vars_, map(float, coefs))))


def _merge(pairs):
    t = {}
    for v, c in pairs:
        t[v] = t.get(v, 0.0) + c
    return t


class Constraint:
    """expr (sense) 0, normalized to  terms (sense) rhs."""

    __slots__ = ("terms", "sense", "rhs", "name", "index")

    def __init__(self, expr, sense):
        self.terms = expr.terms
        self.sense = sense
        self.rhs = -expr.const
        self.name = None
        self.index = None

    def __bool__(self):  # guard against chained comparisons
        raise TypeError("constraints have no truth value")


class Model:
    def __init__(self, name=""):
        self.name = name
        self.sense = "min"  # or "max"
        self.vars: list[Var] = []
        self.constrs: list[Constraint] = []

    def add_var(self, lb=0.0, ub=INF, obj=0.0, vtype="C", name=None):
        vtype = str(vtype)
        if vtype not in ("C", "I", "B"):
            raise ValueError(f"bad vtype {vtype!r}")
        if vtype == "B":
            lb, ub = max(0.0, float(lb)), min(1.0, float(ub))
        v = Var(len(self.vars), name if name is not None else f"C{len(self.vars)}",
                float(lb), float(ub), float(obj), vtype)
        self.vars.append(v)
        return v

    def add_constr(self, con, name=None):
        if not isinstance(con, Constraint):
            raise TypeError("add_constr expects expr <= / >= / == rhs")
        con.index = len(self.constrs)
        con.name = name if name is not None else f"R{con.index}"
        self.constrs.append(con)
        return con

    # ---- statistics
    @property
    def num_vars(self):
        return len(self.vars)

    @property
    def num_constrs(self):
        return len(self.constrs)

    @property
    def num_nonzeros(self):
        return sum(1 for c in self.constrs for a in c.terms.values() if a != 0.0)

    @property
    def num_int_vars(self):
        return sum(1 for v in self.vars if v.vtype != "C")

    # ---- MPS output
    def write(self, path):
        with open(path, "w") as f:
            f.write(self.to_mps())

    def to_mps(self):
        out = [f"NAME {self.name}"]
        if self.sense == "max":
            out.append("OBJSENSE MAX")
        elif self.sense != "min":
            raise ValueError(f"bad sense {self.sense!r}")
        out.append("ROWS")
        out.append(" N  OBJ")
        for c in self.constrs:
            out.append(f" {c.sense}  {c.name:<8}")
        # column-wise matrix (entries in row order)
        cols = [[] for _ in self.vars]
        for c in self.constrs:
            for v, a in c.terms.items():
                if a != 0.0:
                    cols[v.index].append((c.index, a))
        out.append("COLUMNS")
        in_int = False
        for v in self.vars:
            is_int = v.vtype != "C"
            if is_int != in_int:
                out.append("    MARKER    'MARKER'                 "
                           + ("'INTORG'" if is_int else "'INTEND'"))
                in_int = is_int
            entries = sorted(cols[v.index])
            if v.obj != 0.0:
                out.append(f"    {v.name:<8}  {'OBJ':<8}  {_num(v.obj)}")
            elif not entries:
                out.append(f"    {v.name:<8}  {'OBJ':<8}  0")
            for i, a in entries:
                out.append(f"    {v.name:<8}  {self.constrs[i].name:<8}  {_num(a)}")
        if in_int:
            out.append("    MARKER    'MARKER'                 'INTEND'")
        out.append("RHS")
        for c in self.constrs:
            if c.rhs != 0.0:
                out.append(f"    RHS1      {c.name:<8}  {_num(c.rhs)}")
        out.append("BOUNDS")
        for v in self.vars:
            out.extend(_bounds(v))
        out.append("ENDATA")
        return "\n".join(out) + "\n"


def _bound(kind, name, val=None):
    s = f" {kind} BND1      {name:<8}"
    return s if val is None else f"{s}  {_num(val)}"


def _bounds(v):
    lb, ub, n = v.lb, v.ub, v.name
    if v.vtype == "B" and lb == 0.0 and ub == 1.0:
        return [_bound("BV", n)]
    if lb == ub:
        return [_bound("FX", n, lb)]
    out = []
    if v.vtype != "C" and ub == INF:
        # integer with no upper bound: LI keeps ub = +inf under every reader
        out.append(_bound("LI", n, lb) if lb != -INF else _bound("MI", n))
        if lb == -INF:
            out.append(_bound("PL", n))
        return out
    if lb == -INF and ub == INF:
        return [_bound("FR", n)]
    if lb == -INF:
        out.append(_bound("MI", n))
    elif lb != 0.0:
        out.append(_bound("LO", n, lb))
    if ub != INF:
        out.append(_bound("UP", n, ub))
    return out


def _num(x):
    """Shortest exact decimal up to 15 significant digits, else %.16e."""
    x = float(x)
    if x == 0.0:
        return "0"
    s = f"{x:.15g}"
    if float(s) == x:
        return s
    return f"{x:.16e}"
