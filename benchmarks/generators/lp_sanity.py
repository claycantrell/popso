"""Hand-written tiny LP/MIP sanity instances.

These target the LP engine and the MPS reader's edge cases (RANGES, all
BOUNDS types, negative UP semantics, OBJSENSE forms, MARKER blocks, objective
constants, empty rows/columns, extra N rows, fixed format, ...).

Each case is written verbatim as `lpsanity_<case>_0.mps`.  EXPECTED holds
hand-derived answers (status, objective) used by reference.py as an extra
check against the reference solvers; None = not derived by hand.
"""

from __future__ import annotations

from pathlib import Path

CASES: dict[str, str] = {}
EXPECTED: dict[str, tuple] = {}

# Cases whose semantics differ between readers.  We follow the CPLEX
# convention: a negative UP bound on a column without an explicit lower bound
# makes the lower bound -inf.  HiGHS 1.15 and SCIP keep lb = 0 and report
# infeasible (SCIP even fails to read the integer case).  For these the
# reference uses only the hand-derived expectation.
HAND_DERIVED_ONLY = {"negative_ub", "mip_negative_ub"}


def case(name, expected=None):
    def deco(text):
        CASES[name] = text.strip("\n") + "\n"
        EXPECTED[name] = expected
        return text

    return deco


# max 3x + 2y, x + y <= 4, x + 3y <= 6, x <= 3  ->  x=3, y=1, obj 11
case("feasible", ("optimal", 11.0))("""
* simple feasible LP, maximization using an OBJSENSE section
NAME feasible
OBJSENSE
    MAX
ROWS
 N  obj
 L  c1
 L  c2
COLUMNS
    x  obj  3  c1  1
    x  c2  1
    y  obj  2  c1  1
    y  c2  3
RHS
    rhs  c1  4  c2  6
BOUNDS
 UP bnd  x  3
ENDATA
""")

# same model, inline "OBJSENSE MAX" (as written by common modelling tools)
case("objsense_inline", ("optimal", 11.0))("""
NAME objsense_inline
OBJSENSE MAX
ROWS
 N  obj
 L  c1
 L  c2
COLUMNS
    x  obj  3  c1  1
    x  c2  1
    y  obj  2  c1  1
    y  c2  3
RHS
    rhs  c1  4  c2  6
BOUNDS
 UP bnd  x  3
ENDATA
""")

# section form with MAXIMIZE keyword and explicit MIN-like objective name
case("objsense_maximize", ("optimal", 11.0))("""
NAME objsense_maximize
OBJSENSE
    MAXIMIZE
ROWS
 N  profit
 L  c1
 L  c2
COLUMNS
    x  profit  3  c1  1
    x  c2  1
    y  profit  2  c1  1
    y  c2  3
RHS
    rhs  c1  4  c2  6
BOUNDS
 UP bnd  x  3
ENDATA
""")

case("infeasible", ("infeasible", None))("""
* x + y <= 1 and x + y >= 2
NAME infeasible
ROWS
 N  obj
 L  c1
 G  c2
COLUMNS
    x  obj  1  c1  1
    x  c2  1
    y  obj  1  c1  1
    y  c2  1
RHS
    rhs  c1  1  c2  2
ENDATA
""")

case("infeasible_bounds", ("infeasible", None))("""
* x <= 1, y <= 1 (bounds) but x + y >= 3
NAME infeasible_bounds
ROWS
 N  obj
 G  c1
COLUMNS
    x  obj  1  c1  1
    y  obj  1  c1  1
RHS
    rhs  c1  3
BOUNDS
 UP bnd  x  1
 UP bnd  y  1
ENDATA
""")

case("unbounded", ("unbounded", None))("""
* min -x - y, x - y <= 1: ray (1,1)
NAME unbounded
ROWS
 N  obj
 L  c1
COLUMNS
    x  obj  -1  c1  1
    y  obj  -1  c1  -1
RHS
    rhs  c1  1
ENDATA
""")

case("unbounded_free", ("unbounded", None))("""
* min x with x free, x - y <= 3, y in [0, 10]
NAME unbounded_free
ROWS
 N  obj
 L  c1
COLUMNS
    x  obj  1  c1  1
    y  c1  -1
RHS
    rhs  c1  3
BOUNDS
 FR bnd  x
 UP bnd  y  10
ENDATA
""")

# Beale's classic cycling example
case("beale_cycling", ("optimal", -1.25))("""
NAME beale_cycling
ROWS
 N  obj
 L  r1
 L  r2
 L  r3
COLUMNS
    x4  obj  -0.75  r1  0.25
    x4  r2  0.5
    x5  obj  20  r1  -8
    x5  r2  -12
    x6  obj  -0.5  r1  -1
    x6  r2  -0.5  r3  1
    x7  obj  6  r1  9
    x7  r2  3
RHS
    rhs  r3  1
ENDATA
""")

# all six constraints tight at (1,1): primal degenerate vertex
case("degenerate", ("optimal", -2.0))("""
NAME degenerate
ROWS
 N  obj
 L  r1
 L  r2
 L  r3
 L  r4
 L  r5
 L  r6
COLUMNS
    x  obj  -1  r1  1
    x  r3  1  r4  2
    x  r5  1  r6  3
    y  obj  -1  r2  1
    y  r3  1  r4  1
    y  r5  2  r6  1
RHS
    rhs  r1  1  r2  1
    rhs  r3  2  r4  3
    rhs  r5  3  r6  4
ENDATA
""")

# max x + y, x + y <= 4, x <= 3, y <= 3: a whole edge is optimal
case("dual_degenerate", ("optimal", 4.0))("""
NAME dual_degenerate
OBJSENSE MAX
ROWS
 N  obj
 L  c1
COLUMNS
    x  obj  1  c1  1
    y  obj  1  c1  1
RHS
    rhs  c1  4
BOUNDS
 UP bnd  x  3
 UP bnd  y  3
ENDATA
""")

# min 2x + y + z, x - y = 2, x + y - z = 0, x FR, y MI+PL (free), z >= 0
# -> x = 1, y = -1, z = 0, obj 1
case("freevars", ("optimal", 1.0))("""
NAME freevars
ROWS
 N  obj
 E  e1
 E  e2
COLUMNS
    x  obj  2  e1  1
    x  e2  1
    y  obj  1  e1  -1
    y  e2  1
    z  obj  1  e2  -1
RHS
    rhs  e1  2
BOUNDS
 FR bnd  x
 MI bnd  y
 PL bnd  y
ENDATA
""")

# Ranged rows of every flavour:
#   r1: L, rhs 10, R  4  ->  6 <= x + y <= 10
#   r2: G, rhs  2, R  3  ->  2 <= x - y <= 5
#   r3: E, rhs  3, R  2  ->  3 <= y + z <= 5
#   r4: E, rhs  4, R -3  ->  1 <= z <= 4
#   r5: L, rhs 12, R -5  ->  7 <= x + z <= 12   (sign of R ignored for L/G)
# min x + 2y + 3z -> 13 (e.g. x = 6, y = 2, z = 1; a whole edge is optimal)
case("ranges", ("optimal", 13.0))("""
NAME ranges
ROWS
 N  obj
 L  r1
 G  r2
 E  r3
 E  r4
 L  r5
COLUMNS
    x  obj  1  r1  1
    x  r2  1  r5  1
    y  obj  2  r1  1
    y  r2  -1  r3  1
    z  obj  3  r3  1
    z  r4  1  r5  1
RHS
    rhs  r1  10  r2  2
    rhs  r3  3  r4  4
    rhs  r5  12
RANGES
    rng  r1  4  r2  3
    rng  r3  2  r4  -3
    rng  r5  -5
ENDATA
""")

# Negative UP on a column with default lower bound: lb becomes -inf.
# min x + y, x + y >= -7, x <= -2, y <= -1   -> obj -7
# (a reader keeping lb = 0 would report infeasible)
case("negative_ub", ("optimal", -7.0))("""
NAME negative_ub
ROWS
 N  obj
 G  c1
COLUMNS
    x  obj  1  c1  1
    y  obj  1  c1  1
RHS
    rhs  c1  -7
BOUNDS
 UP bnd  x  -2
 UP bnd  y  -1
ENDATA
""")

# Negative lower bounds: min x + 2y, x in [-3, 5], y >= -10,
# x + y >= -8, x - y <= 4 -> y >= x - 4 ; best: y as small as possible:
# y >= -8 - x and y >= x - 4 and y >= -10. x=-2 gives y >= -6 both -> obj -14;
# objective x + 2y with y = max(-8-x, x-4): for x <= -2 y=-8-x -> obj -16 - x, min at x=-2 -> -14
# for x >= -2 y = x - 4 -> obj 3x - 8, min at x=-2 -> -14.
case("negative_lb", ("optimal", -14.0))("""
NAME negative_lb
ROWS
 N  obj
 G  c1
 L  c2
COLUMNS
    x  obj  1  c1  1
    x  c2  1
    y  obj  2  c1  1
    y  c2  -1
RHS
    rhs  c1  -8  c2  4
BOUNDS
 LO bnd  x  -3
 UP bnd  x  5
 LO bnd  y  -10
ENDATA
""")

# Variables with only an upper bound.
# x: UP 5 (lb stays 0); y: MI then UP 7 -> (-inf, 7]; z: MI only -> (-inf, +inf)
# max x + y - z, z >= -4 (row), x + y <= 100 -> 5 + 7 + 4 = 16
case("upper_only", ("optimal", 16.0))("""
NAME upper_only
OBJSENSE
    MAX
ROWS
 N  obj
 L  c1
 G  c2
COLUMNS
    x  obj  1  c1  1
    y  obj  1  c1  1
    z  obj  -1  c2  1
RHS
    rhs  c1  100  c2  -4
BOUNDS
 UP bnd  x  5
 MI bnd  y
 UP bnd  y  7
 MI bnd  z
ENDATA
""")

# FX bounds, including a negative fixed value.
# x = 2.5, y = -1, w >= x - y = 3.5; min x + y + w = 5
case("fixed", ("optimal", 5.0))("""
NAME fixed
ROWS
 N  obj
 G  c1
COLUMNS
    x  obj  1  c1  -1
    y  obj  1  c1  1
    w  obj  1  c1  1
RHS
    rhs  c1  0
BOUNDS
 FX bnd  x  2.5
 FX bnd  y  -1
ENDATA
""")

# RHS on the objective row is the NEGATED objective constant.
# min x + 10, x >= 3 -> 13
case("objoffset", ("optimal", 13.0))("""
NAME objoffset
ROWS
 N  obj
 G  c1
COLUMNS
    x  obj  1  c1  1
RHS
    rhs  obj  -10
    rhs  c1  3
ENDATA
""")

# Column x appears only in the objective, z has no entries except a zero
# objective coefficient.  min -x + y, x <= 4, y >= 1 -> -3
case("empty_col", ("optimal", -3.0))("""
NAME empty_col
ROWS
 N  obj
 G  c1
COLUMNS
    x  obj  -1
    y  obj  1  c1  1
    z  obj  0
RHS
    rhs  c1  1
BOUNDS
 UP bnd  x  4
 FR bnd  z
ENDATA
""")

# Row with no coefficients (0 <= 5) is harmless.  min x + y, x + y >= 2 -> 2
case("empty_row", ("optimal", 2.0))("""
NAME empty_row
ROWS
 N  obj
 L  empty
 G  c1
COLUMNS
    x  obj  1  c1  1
    y  obj  1  c1  1
RHS
    rhs  empty  5  c1  2
ENDATA
""")

# Row with no coefficients and 0 >= 1: infeasible.
case("empty_row_infeasible", ("infeasible", None))("""
NAME empty_row_infeasible
ROWS
 N  obj
 G  empty
 G  c1
COLUMNS
    x  obj  1  c1  1
RHS
    rhs  empty  1  c1  2
ENDATA
""")

# A second N row must be ignored (it is not a constraint).
# min x + y, x + y >= 2, x - y = 0 -> 2
case("extra_free_row", ("optimal", 2.0))("""
NAME extra_free_row
ROWS
 N  obj
 N  cost2
 G  c1
 E  c2
COLUMNS
    x  obj  1  cost2  -100
    x  c1  1  c2  1
    y  obj  1  cost2  50
    y  c1  1  c2  -1
RHS
    rhs  c1  2
ENDATA
""")

# Linearly dependent equalities (rank deficient).
# x + y + z = 3, 2x + 2y + 2z = 6, x - y = 0, min x + 2y + 3z = 9 - 3x -> x = 1.5 -> 4.5
case("redundant_eq", ("optimal", 4.5))("""
NAME redundant_eq
ROWS
 N  obj
 E  e1
 E  e2
 E  e3
COLUMNS
    x  obj  1  e1  1
    x  e2  2  e3  1
    y  obj  2  e1  1
    y  e2  2  e3  -1
    z  obj  3  e1  1
    z  e2  2
RHS
    rhs  e1  3  e2  6
ENDATA
""")

case("inconsistent_eq", ("infeasible", None))("""
* x + y = 3 and 2x + 2y = 7
NAME inconsistent_eq
ROWS
 N  obj
 E  e1
 E  e2
COLUMNS
    x  obj  1  e1  1
    x  e2  2
    y  obj  1  e1  1
    y  e2  2
RHS
    rhs  e1  3  e2  7
ENDATA
""")

# Klee-Minty cube, n = 6: max sum 2^(n-j) x_j,
#   sum_{j<i} 2^(i-j+1) x_j + x_i <= 5^i  ->  optimum 5^6 = 15625
def _klee_minty(n=6):
    rows = [f" L  r{i}" for i in range(1, n + 1)]
    cols = []
    for j in range(1, n + 1):
        cols.append(f"    x{j}  obj  {2 ** (n - j)}")
        for i in range(j, n + 1):
            coef = 1 if i == j else 2 ** (i - j + 1)
            cols.append(f"    x{j}  r{i}  {coef}")
    rhs = [f"    rhs  r{i}  {5 ** i}" for i in range(1, n + 1)]
    return "\n".join(
        ["NAME klee_minty", "OBJSENSE", "    MAX", "ROWS", " N  obj", *rows,
         "COLUMNS", *cols, "RHS", *rhs, "ENDATA"]
    )


case("klee_minty", ("optimal", 15625.0))(_klee_minty())

# Badly scaled coefficients (1e-4 .. 1e5).
case("illscaled", ("optimal", None))("""
NAME illscaled
ROWS
 N  obj
 G  c1
 L  c2
 G  c3
COLUMNS
    x  obj  1.0E+00  c1  1.0E-04
    x  c2  2.5E+00  c3  1
    y  obj  3.0E+03  c1  1.0E+01
    y  c2  1.0E+05
    z  obj  2.0E-03  c1  1.0E+00
    z  c3  1.0E+04
RHS
    rhs  c1  1.5E+00  c2  2.0E+05
    rhs  c3  2.5E+04
BOUNDS
 UP bnd  z  1.0E+02
ENDATA
""")

# Fixed-format MPS (fields in columns 2-3, 5-12, 15-22, 25-36, 40-47, 50-61).
case("fixed_format", ("optimal", None))("""
NAME          FIXEDFMT
ROWS
 N  COST
 L  LIM1
 G  LIM2
 E  MYEQN
COLUMNS
    XONE      COST         1.0   LIM1         1.0
    XONE      LIM2         1.0
    YTWO      COST         2.0   LIM1         1.0
    YTWO      MYEQN       -1.0
    ZTHREE    COST        -1.0   MYEQN        1.0
RHS
    RHS       LIM1         4.0   LIM2         1.0
    RHS       MYEQN        7.0
BOUNDS
 UP BND       XONE         4.0
 LO BND       YTWO        -1.0
 UP BND       YTWO         1.0
ENDATA
""")

# ---------------------------------------------------------------- MIP cases

# MARKER block integers with explicit bounds, BV, and LI/UI outside MARKER.
#   a: BV; b: MARKER, UP 3; c: LI -2 / UI 2 (integer via bound type);
#   d: continuous <= 1.5
# max 5a + 4b + 3c + d  s.t.  2a + 3b + c + d <= 7.5,  a + b - c >= 1
case("mip_marker", ("optimal", None))("""
NAME mip_marker
OBJSENSE MAX
ROWS
 N  obj
 L  c1
 G  c2
COLUMNS
    MARKER  'MARKER'  'INTORG'
    a  obj  5  c1  2
    a  c2  1
    b  obj  4  c1  3
    b  c2  1
    MARKER  'MARKER'  'INTEND'
    c  obj  3  c1  1
    c  c2  -1
    d  obj  1  c1  1
RHS
    rhs  c1  7.5  c2  1
BOUNDS
 BV bnd  a
 UP bnd  b  3
 LI bnd  c  -2
 UI bnd  c  2
 UP bnd  d  1.5
ENDATA
""")

# Integer column with negative UP: lb becomes -inf.
# min x, x >= -10.5 (row), x integer <= -3 -> -10
case("mip_negative_ub", ("optimal", -10.0))("""
NAME mip_negative_ub
ROWS
 N  obj
 G  c1
COLUMNS
    MARKER  'MARKER'  'INTORG'
    x  obj  1  c1  1
    MARKER  'MARKER'  'INTEND'
RHS
    rhs  c1  -10.5
BOUNDS
 UP bnd  x  -3
ENDATA
""")

# LP relaxation feasible, but 2x + 2y = 3 has no integer solution.
case("mip_parity_infeasible", ("infeasible", None))("""
NAME mip_parity_infeasible
ROWS
 N  obj
 E  c1
COLUMNS
    MARKER  'MARKER'  'INTORG'
    x  obj  1  c1  2
    y  obj  1  c1  2
    MARKER  'MARKER'  'INTEND'
RHS
    rhs  c1  3
BOUNDS
 UP bnd  x  10
 UP bnd  y  10
ENDATA
""")

# max x + y, |x - y| <= 2.5, x, y integer >= 0: unbounded along (1,1).
# The PL bounds matter: without any BOUNDS entry MARKER integers are binary.
case("mip_unbounded", ("unbounded", None))("""
NAME mip_unbounded
OBJSENSE MAX
ROWS
 N  obj
 L  c1
 L  c2
COLUMNS
    MARKER  'MARKER'  'INTORG'
    x  obj  1  c1  1
    x  c2  -1
    y  obj  1  c1  -1
    y  c2  1
    MARKER  'MARKER'  'INTEND'
RHS
    rhs  c1  2.5  c2  2.5
BOUNDS
 PL bnd  x
 PL bnd  y
ENDATA
""")

# MARKER-integer default bounds (convention shared by HiGHS/SCIP and others):
#   x: no BOUNDS entry   -> binary [0, 1]
#   y: only LO 2         -> [2, +inf]
#   z: only MI           -> (-inf, +inf)
# max x + y - z, x + y <= 7.5, z >= -3.5 -> x + y = 7, z = -3 -> 10
case("mip_marker_defaults", ("optimal", 10.0))("""
NAME mip_marker_defaults
OBJSENSE MAX
ROWS
 N  obj
 L  c1
 G  c2
COLUMNS
    MARKER  'MARKER'  'INTORG'
    x  obj  1  c1  1
    y  obj  1  c1  1
    z  obj  -1  c2  1
    MARKER  'MARKER'  'INTEND'
RHS
    rhs  c1  7.5  c2  -3.5
BOUNDS
 LO bnd  y  2
 MI bnd  z
ENDATA
""")

# Equality knapsack with general integers:
# 3x + 5y + 7z = 101, min x + y + z -> z=13 (91), 10 left = 5*2 -> 15
case("mip_equality", ("optimal", 15.0))("""
NAME mip_equality
ROWS
 N  obj
 E  c1
COLUMNS
    MARKER  'MARKER'  'INTORG'
    x  obj  1  c1  3
    y  obj  1  c1  5
    z  obj  1  c1  7
    MARKER  'MARKER'  'INTEND'
RHS
    rhs  c1  101
BOUNDS
 UP bnd  x  100
 UP bnd  y  100
 UP bnd  z  100
ENDATA
""")

# Small 0-1 knapsack with BV bounds only (no MARKER).
# max 10a + 13b + 7c + 8d + 4e, 4a + 6b + 3c + 4d + 2e <= 11
case("mip_knapsack", ("optimal", None))("""
NAME mip_knapsack
OBJSENSE MAX
ROWS
 N  obj
 L  cap
COLUMNS
    a  obj  10  cap  4
    b  obj  13  cap  6
    c  obj  7  cap  3
    d  obj  8  cap  4
    e  obj  4  cap  2
RHS
    rhs  cap  11
BOUNDS
 BV bnd  a
 BV bnd  b
 BV bnd  c
 BV bnd  d
 BV bnd  e
ENDATA
""")

# Mixed: free integer, fixed integer, continuous; objective constant.
# min x - y + 2w + 100 (RHS obj -100)
#   x integer free, y integer in [-5, 5], w continuous >= 0, v integer fixed 3
#   x + y >= -3.5,  x - y + w >= 0.5,  x <= 4 (row), v + x <= 6
case("mip_mixed_bounds", ("optimal", None))("""
NAME mip_mixed_bounds
ROWS
 N  obj
 G  c1
 G  c2
 L  c3
 L  c4
COLUMNS
    MARKER  'MARKER'  'INTORG'
    x  obj  1  c1  1
    x  c2  1  c3  1
    x  c4  1
    y  obj  -1  c1  1
    y  c2  -1
    v  c4  1
    MARKER  'MARKER'  'INTEND'
    w  obj  2  c2  1
RHS
    rhs  obj  -100
    rhs  c1  -3.5  c2  0.5
    rhs  c3  4  c4  6
BOUNDS
 FR bnd  x
 LO bnd  y  -5
 UP bnd  y  5
 FX bnd  v  3
ENDATA
""")


def write_all(outdir: Path) -> list[Path]:
    outdir.mkdir(parents=True, exist_ok=True)
    paths = []
    for name, text in CASES.items():
        p = outdir / f"lpsanity_{name}_0.mps"
        p.write_text(text)
        paths.append(p)
    return paths


def hand_derived_only_instances() -> set:
    return {f"lpsanity_{k}_0" for k in HAND_DERIVED_ONLY}


def expected_by_instance() -> dict:
    return {f"lpsanity_{k}_0": v for k, v in EXPECTED.items() if v is not None}
