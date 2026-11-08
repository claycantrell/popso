# popso

**popso** is an open-source solver for linear and mixed-integer linear programs (LP / MILP), convex quadratic
programs (QP / MIQP), and convex quadratically constrained and second-order cone programs (QCQP / MIQCQP,
SOCP / MISOCP). It is written from scratch in C++20.

`pip install popso` gives you both:

- the **`popso` command**: `popso model.mps --time-limit 60` (MPS, LP and CBF files, plain or gzipped);
- the **Python API**:

```python
import popso

m = popso.Model("example")
x = m.add_var(ub=10, vtype=popso.INTEGER, name="x")
y = m.add_var(name="y")
m.add_constr(x + 2 * y <= 14)
m.add_constr(3 * x - y >= 0)
m.set_objective(3 * x + 4 * y, sense=popso.MAXIMIZE)
r = m.solve(time_limit=10)
print(r.status, r.objective, r[x], r[y])
```

Documentation, benchmarks against HiGHS and SCIP, and the changelog are on
[GitHub](https://github.com/gozwei/popso). License: Apache-2.0.
