"""popso: Python bindings (ctypes) for the popso LP/MILP/QP/MIQP solver.

    import popso
    m = popso.Model("example")
    x = m.add_var(ub=10, vtype=popso.INTEGER, name="x")
    y = m.add_var(name="y")
    m.add_constr(x + 2 * y <= 14)
    m.add_constr(3 * x - y >= 0)
    m.set_objective(3 * x + 4 * y, sense=popso.MAXIMIZE)
    r = m.solve(time_limit=10)
    print(r.status, r.objective, r[x], r["y"])

The shared library is searched in $POPSO_LIBRARY (a file or a directory), next to this
package, in the build/ directory of a source checkout, and on the system library path.
"""

from ._lib import path as library_path
from .model import (
    AT_LOWER,
    AT_UPPER,
    AT_ZERO,
    BASIC,
    BINARY,
    CONTINUOUS,
    INF,
    INTEGER,
    MAXIMIZE,
    MINIMIZE,
    PARAMS,
    SEMICONTINUOUS,
    SEMIINTEGER,
    SOS1,
    SOS2,
    Constr,
    Constraint,
    LinExpr,
    Model,
    Params,
    PopsoError,
    QuadExpr,
    Result,
    Var,
    read_lp,
    read_mps,
    version,
)

__version__ = version()

__all__ = [
    "AT_LOWER", "AT_UPPER", "AT_ZERO", "BASIC", "BINARY", "CONTINUOUS", "INF", "INTEGER", "MAXIMIZE", "MINIMIZE", "PARAMS", "SEMICONTINUOUS", "SEMIINTEGER",
    "SOS1", "SOS2", "Constr", "Constraint", "LinExpr", "Model", "Params", "PopsoError", "QuadExpr", "Result", "Var",
    "library_path", "read_lp", "read_mps", "version",
]
