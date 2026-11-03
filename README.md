<p align="center"><img src="docs/popso-banner.png" alt="popso" width="640"></p>

# popso

[![ci](https://github.com/gozwei/popso/actions/workflows/ci.yml/badge.svg)](https://github.com/gozwei/popso/actions/workflows/ci.yml)

**popso** (Polkowski's Optimization Solver) is an open-source solver for **mixed-integer linear programs (MILP)** and **linear programs
(LP)**, and, since v0.4, for **convex quadratic programs (QP)** and **mixed-integer convex quadratic programs (MIQP)**. It is written from scratch in C++20 using only the standard library.

This is a fun side project by [Marcin Polkowski](https://marcinpolkowski.com). The goal is
a fully open-source solver that holds its own against the best solvers on real
problems. It is built with heavy use of AI coding agents. More problem classes are planned.

![MIPLIB: instances solved vs. time](docs/miplib_core.svg)

On a 63-instance [MIPLIB 2017](https://miplib.zib.de) subset it currently solves more
instances, and is faster overall, than the open-source solvers HiGHS and SCIP:

| solver | solved | shifted geomean (s) |
|---|---|---|
| popso | 59/63 | 10.5 |
| HiGHS 1.15 | 53/63 | 18.9 |
| SCIP 10 | 50/63 | 21.3 |

On the 43 harder instances of the same collection that no solver finishes in 2 minutes, with a
10-minute limit, popso solves 9, HiGHS 12 and SCIP 6.

These are single-machine numbers with caveats. See **[BENCHMARKS.md](BENCHMARKS.md)** for the
setup, per-instance results, the hard instances, the generated benchmark families and popso's
progress across versions, and
**[CHANGELOG.md](CHANGELOG.md)** for what each version can do.

## Build

You need CMake ≥ 3.20 and a C++20 compiler (clang, gcc or MSVC). No other dependencies.

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DPOPSO_BUILD_TESTS=ON
cmake --build build -j 4
./build/popso_tests            # unit tests
```

On Windows with Visual Studio, add `--config Release` to the build command; the binaries are
then in `build/Release/`.

The solver is developed and benchmarked on macOS (Apple silicon). Continuous integration
builds it and runs the unit tests and a small end-to-end check on Linux (gcc and clang),
macOS (clang) and Windows (MSVC).

## Usage

```sh
./build/popso model.mps --sol model.sol --time-limit 60 --threads 4
```

The solver reads a model in MPS format (free or fixed) or in the CPLEX LP format and prints
a summary. Besides continuous, integer and binary variables, models may have
semi-continuous and semi-integer variables, SOS1/SOS2 constraints and indicator
constraints (MPS `SC`/`SI` bounds, `SOS` and `INDICATORS` sections; LP `Semi-Continuous` and
`SOS` sections and `b = 1 -> ...` constraints). The objective may have a convex quadratic part
(see [Quadratic objectives](#quadratic-objectives-qp-and-miqp) below). The format is chosen by the file extension (`.lp` is LP, anything else MPS,
case-insensitive); `--format lp|mps` overrides it:

```
Status: OPTIMAL
Objective: 12574
Bound: 12574
Gap: 0.0000%
Nodes: 1
Time: 0.003s
```

With `--sol` it also writes the solution as a plain `name value` file:

```
# Solution for model knapsack_small_1
# Status = OPTIMAL
# Objective value = 12574
x0 0
x1 1
...
```

| option | meaning |
|---|---|
| `--sol <file>` | write the solution file |
| `--start <file>` | starting solution in the same `name value` format; may be partial (see below) |
| `--duals <file>` | write row duals and reduced costs (`row <name> <value>` / `col <name> <value>` lines; continuous QPs) |
| `--format <mps\|lp>` | input format (default: from the file extension) |
| `--time-limit <sec>` | wall-clock time limit (default: none) |
| `--threads <n>` | number of threads (0 = all cores, the default) |
| `--gap <rel>` | relative optimality gap at which to stop (default 1e-4) |
| `--deterministic` | reproducible parallel search (slower; the default is opportunistic) |
| `--seed <n>` | random seed |
| `--racing <k>` | run k differently seeded searches at the same time |
| `--mem-limit <MB>` | approximate memory cap for the search tree (default 4096) |
| `--node-limit <n>` | stop after n branch-and-bound nodes |
| `--no-presolve` | disable presolve |
| `--lp-method <m>` | LP algorithm for a pure LP and the MIP root LP: `auto` (default), `simplex` or `barrier` |
| `--verbosity <0-3>` | 0 silent, 1 summary (default), 2 progress, 3 debug |
| `--version` | print the version |

**Starting solutions.** `--start` takes a solution file, for example one written by an earlier run
with `--sol`. A complete feasible start becomes the first incumbent. If only some values are given, or
the start is infeasible, the given integer values are fixed and the rest is completed by a short
sub-MIP. If that fails, the start is ignored and the log says so.

### Quadratic objectives (QP and MIQP)

The objective can be `c'x + 1/2 x'Qx` with a symmetric Q. In MPS, Q comes in a `QUADOBJ` section
(one triangle; an off-diagonal entry stands for both Q_ij and Q_ji) or in a `QMATRIX` / CPLEX
`QSECTION` section (the whole matrix). In the LP format, the quadratic terms go in brackets in the
objective:

```
Minimize
 obj: - x - 2 y + [ 2 x ^ 2 + 2 x * y + 4 y ^ 2 ] / 2
Subject To
 c1: x + y <= 3
End
```

- The objective must be **convex**: Q positive semidefinite when minimizing, negative
  semidefinite when maximizing. popso checks this before solving (a sparse factorization of the
  scaled Q) and reports the status `NONCONVEX` otherwise; Q that is semidefinite only up to a
  relative 1e-6 (round-off in the input) is accepted with that much added to its diagonal.
- **Presolve** for quadratic models removes fixed and empty columns (their quadratic terms move
  into the linear part), fixes columns without quadratic terms by dual arguments, turns
  singleton rows into bounds and drops empty and redundant rows; no column is substituted.
- **Continuous QPs** are solved by the barrier method with Q in the Newton system (the augmented
  system; the normal equations when Q is diagonal). There is no crossover and no basis: the
  result is the primal point, the objective and the **dual values** (row duals y and reduced
  costs d = c + Qx - A'y of the final interior point, so complementarity holds to the barrier's
  tolerance; when minimizing, y >= 0 at an active lower row bound and <= 0 at an upper one, the
  opposite when maximizing). Infeasible and unbounded QPs are confirmed by the simplex method on
  auxiliary LPs.
- **MIQPs** (integer variables, linear constraints, convex objective) are solved by
  branch-and-bound on QP relaxations: the root by the barrier, the nodes by a dual active-set
  method warm-started from the parent node. Node bounds are Lagrangian bounds that stay valid for
  inexact solves, and they also tighten integer bounds (reduced-cost fixing). The search uses
  domain propagation, reliability branching, rounding, a QP dive and feasibility jump on a
  linearized objective; it is single-threaded and deterministic.
- Cutting planes are not used for quadratic models yet, and quadratic constraints are not
  supported (MPS `QCMATRIX`, LP brackets in constraints: an error).

The C API sets Q with `popso_set_q` / `popso_add_q_entries` (and reads it with
`popso_get_q`) and returns the duals of a continuous QP with `popso_get_duals` /
`popso_get_reduced_costs` (`popso_has_duals`); in Python, products of variables in
`set_objective` make a quadratic objective, and `Result.duals` / `Result.reduced_costs` (or
`r.dual(constr)`, `r.reduced_cost(var)`) hold the duals:

```python
x = m.add_var(lb=-10, ub=10, name="x")
y = m.add_var(lb=-10, ub=10, vtype=popso.INTEGER, name="y")
m.set_objective((x - 1) ** 2 + x * y + 2 * y * y - 3 * y)   # a popso.QuadExpr
r = m.solve()
```

### C library

The build also produces a shared library (`build/libpopso.dylib`, `libpopso.so` or `popso.dll`) with
a C API declared in [`include/popso.h`](include/popso.h):

```c
#include "popso.h"

popso_model* m = popso_create();
int x, y;
popso_add_var(m, 0.0, 10.0, 3.0, POPSO_INTEGER, "x", &x);
popso_add_var(m, 0.0, POPSO_INFINITY, 4.0, POPSO_CONTINUOUS, "y", &y);
int cols[2] = {x, y};
double vals[2] = {1.0, 2.0};
popso_add_row(m, 2, cols, vals, -POPSO_INFINITY, 14.0, "c1", NULL);  /* x + 2y <= 14 */
popso_set_objective_sense(m, POPSO_MAXIMIZE);
popso_set_param(m, "time_limit", 10.0);
if (popso_solve(m) == POPSO_OK && popso_has_solution(m)) {
  double obj, sol[2];
  popso_get_objective(m, &obj);
  popso_get_solution(m, sol, 2);
}
popso_free(m);
```

Models can also be read with `popso_read_mps` or `popso_read_lp`. `popso_add_sos` and
`popso_add_indicator` add SOS1/SOS2 and indicator constraints, and the variable types
`POPSO_SEMICONTINUOUS` and `POPSO_SEMIINTEGER` make a variable 0 or within its bounds. Functions return error codes, and
`popso_last_error` describes the last failure. Separate models can be solved from different threads
at the same time. The library prints nothing unless you set the `verbosity` parameter.

### Python

The Python package in [`python/`](python) uses `ctypes` to load the shared library. It needs no
compiled extension and no third-party packages.

```python
import popso

m = popso.Model("example")
x = m.add_var(ub=10, vtype=popso.INTEGER, name="x")
y = m.add_var(name="y")
m.add_constr(x + 2 * y <= 14)
m.add_constr(3 * x - y >= 0)
m.set_objective(3 * x + 4 * y, sense=popso.MAXIMIZE)
m.set_start({x: 4})             # optional, may be partial
r = m.solve(time_limit=10, threads=4)
print(r.status, r.objective, r.bound, r[x], r["y"], r.nodes, r.time)

m2 = popso.read_mps("model.mps")

# Logical constraints
b = m.add_var(vtype=popso.BINARY, name="b")
z = m.add_var(lb=2, ub=8, vtype=popso.SEMICONTINUOUS, name="z")   # z = 0 or 2 <= z <= 8
m.add_indicator(b, 0, z <= 0)                                     # b = 0 -> z <= 0
m.add_sos([x, y, z], weights=[1, 2, 3], type=popso.SOS1)          # at most one nonzero
```

To use it from a source checkout, add `python/` to `PYTHONPATH`. The library is then found in
`build/`; you can also point `POPSO_LIBRARY` at it. To install the package, copy the library into the
package and run pip:

```sh
cp build/libpopso.dylib python/popso/     # libpopso.so on Linux, popso.dll on Windows
pip install ./python
```

## What's inside

- **QP and MIQP:** convexity test by a sparse factorization of the scaled Hessian; the barrier
  method with the Hessian in the augmented system; branch-and-bound with a warm-started dual
  active-set QP solver at the nodes, safe Lagrangian bounds, reliability branching and QP-based
  heuristics.
- **LP:** bounded dual simplex with hyper-sparse LU factorization, Forrest–Tomlin updates and
  dual steepest-edge pricing, plus a primal simplex for cleanup. An interior-point (barrier)
  method: Mehrotra predictor-corrector with Gondzio correctors, AMD ordering and a supernodal
  sparse Cholesky factorization (normal equations or augmented system), and a crossover to an
  optimal basis. Large pure LPs run both methods at the same time by default.
- **Logical constraints:** semi-continuous and semi-integer variables, SOS1/SOS2 and
  indicator constraints, reformulated with binaries and big-M rows whose coefficients come
  from bound propagation; without finite bounds, branching on the binary enforces the
  constraint in the subtree.
- **Presolve:** singleton rows and columns, bound, coefficient and GCD tightening, parallel rows
  and columns, dual fixing and dual bound strengthening, dominated columns, aggregation and
  substitution, probing, clique merging and clique extension, with full postsolve.
- **Cutting planes:**
  - Gomory mixed-integer cuts;
  - c-MIR cuts on aggregated rows;
  - lifted knapsack covers;
  - clique, implied-bound and zero-half cuts;
  - tableau cuts at tree nodes;
  - a dynamic cut pool.
- **Primal heuristics:**
  - simple rounding, diving and local search;
  - the feasibility pump and feasibility jump;
  - large-neighbourhood sub-MIPs: RINS, RENS and crossover.
- **Branch-and-bound:**
  - reliability branching on pseudocosts with strong branching;
  - node-level domain propagation and conflict analysis;
  - reduced-cost fixing;
  - symmetry handling (orbital branching, orbitopes, symmetry-breaking reductions from a
    stabilizer chain);
  - root and in-tree restarts.
- **Parallelism:** multithreaded tree search, either opportunistic or fully deterministic
  (`--deterministic` gives the same result on every run with the same thread count).
- **Robustness:** every result is checked against reference solvers. The search tree's
  memory use is bounded. A numerical failure in a node LP is recovered from instead of
  dropping that part of the tree.

## Current limitations

- **Problem classes:** LP, MILP, and convex QP / MIQP (quadratic objective, linear
  constraints). No quadratic constraints, no general nonlinear constraints, no nonconvex
  objectives.
- **Quadratic models:** a basic presolve only, no cutting planes, no parallel MIQP search; the
  node QP solver keeps a dense copy of the Hessian, so MIQPs with more than 3000 columns fall back
  to the much slower barrier at every node. Dual
  values only for continuous QPs (interior point, no basis). MIQP performance is far from the commercial solvers.
- **Model features:** semi-continuous variables, SOS constraints and indicator constraints
  are reformulated with binary variables (big-M rows where bounds allow). There is no
  specialized branching on SOS sets, and presolve is skipped for models with constraints
  that have no finite big-M. No other general constraints (abs, min/max, piecewise-linear
  objectives).
- **Input and interfaces:**
  - MPS and LP-format input only (no gzipped files).
  - The C API and the Python bindings have no callbacks. They return no dual values or reduced
    costs for LPs and MILPs (only for continuous QPs), and a running solve cannot be interrupted.
- **LP algorithms:** the barrier solver is new and not tuned yet; its factorization is slower
  than in mature solvers. It has no homogeneous self-dual model: when it only suspects an
  infeasible or unbounded LP, the simplex confirms it. MIP root LPs use the dual simplex unless
  `--lp-method barrier` is given.
- **Maturity:** this is young software. It is tuned on a limited benchmark set, and
  performance varies from instance to instance and with the random seed. Expect rough edges.

## Repository layout

```
src/milp/      the solver: core/ io/ lp/ presolve/ cuts/ heuristics/ search/ solver/
include/       public C API header (popso.h)
src/api/       C API implementation (shared library popso)
python/        Python bindings (ctypes) and their tests
apps/          command-line front end
tests/         unit tests (own small harness, standard library only) and C API tests
benchmarks/    instance generators, reference results, benchmark runner (Python)
docs/          plots
```

## License

Apache License 2.0 — see [LICENSE](LICENSE). Copyright 2026 Marcin Polkowski.

## Contact

Questions, ideas or interesting models to try? Get in touch via
[marcinpolkowski.com](https://marcinpolkowski.com).
