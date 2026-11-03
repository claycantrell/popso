# Changelog

## Unreleased

## v0.4.0 — 2026-11

New problem classes: convex quadratic programs (QP) and mixed-integer convex quadratic programs
(MIQP): a quadratic objective `c'x + 1/2 x'Qx` with linear constraints. Also better c-MIR cuts
and symmetry-breaking reductions for MILP.

**Benchmarks**
- 63-instance MIPLIB 2017 subset (120 s, 4 threads): 59 solved, shifted geometric mean 10.5 s
  (v0.3.0: 60, 11.4 s); over three seeds 10.5 / 10.9 / 11.4 s against 11.4 / 11.4 / 11.1 s.
  The instance count differs by two instances just above the time limit.
- Generated models, now 148 including 18 QPs and MIQPs (1 thread): 147 solved; SCIP solves 134,
  HiGHS 125 of the 130 it handles (no SOS, indicator or MIQP).
- Maros–Meszaros convex QP test set: all 97 convex problems of the 98 smallest, 32 of the 40
  larger ones (60 s).

**Quadratic objectives**
- Model: a symmetric Q stored as its lower triangle. Input: MPS `QUADOBJ` (one triangle),
  `QMATRIX` and CPLEX `QSECTION` (whole matrix); LP-format objective terms
  `[ 2 x ^ 2 + 3 x * y ] / 2`. The MPS writer writes `QUADOBJ`. Quadratic constraints (MPS
  `QCMATRIX` or a `QSECTION` of another row, brackets in LP constraints) are rejected with a
  "not supported yet" error.
- Convexity check before solving (sparse Cholesky of the unit-diagonal scaled Q plus a tiny
  shift; a negative pivot means a negative eigenvalue). A nonconvex objective (concave
  objective when maximizing) gets the new status `NONCONVEX` (C API `POPSO_STATUS_NONCONVEX`).
- Continuous QPs: the barrier method with Q in the Newton system (augmented system, or the normal
  equations when Q is diagonal), no crossover; infeasibility and unboundedness are confirmed
  by the simplex on auxiliary LPs. Once the dual regularization's error dominates the primal
  residual, the Newton systems are refined against the unregularized matrix (problems with
  large optimal duals such as YAO stalled otherwise).
- Dual values of continuous QPs: row duals and reduced costs (d = c + Qx - A'y) of the
  barrier's final point, in the model's objective sense. CLI `--duals <file>` (`row <name>
  <value>` / `col <name> <value>` lines), C API `popso_has_duals`, `popso_get_duals`,
  `popso_get_reduced_costs`, Python `Result.duals`, `Result.reduced_costs`, `Result.dual(c)`,
  `Result.reduced_cost(x)`. Not for LPs, MILPs and MIQPs.
- Presolve for quadratic models (QP and MIQP): fixed and empty columns (their quadratic terms
  folded into the costs and the objective constant), dual fixing of columns without quadratic
  terms, empty, free and singleton rows, redundant row sides and infeasibility from row
  activities, integer bounds from activities (MIQP). No column substitution, so the Hessian
  stays as it is; the duals are carried back (a singleton row whose bound is active takes the
  reduced cost of its column). The convexity test runs on the presolved Hessian (a cross term
  with a fixed column is linear).
- MIQPs: branch-and-bound (best bound with plunging, domain propagation of the rows),
  deterministic. The root relaxation is solved by the barrier; the nodes by a dual active-set
  method (Goldfarb–Idnani on the KKT system of the working set: dense Cholesky of the free
  columns' Hessian plus a Schur complement of the active rows, later working-set changes as
  borders of that factorization), warm-started from the parent's working set and, when it is
  still cached, its factorization; a singular Hessian gets a tiny proximal term. Node bounds
  are safe Lagrangian bounds (valid for any primal point and row multipliers by convexity),
  which also give reduced-cost tightening of integer bounds at the nodes and at the root.
  Reliability branching (strong branching by the active-set method on candidates with fewer
  than 4 pseudocost observations), plunging into the child with the smaller pseudocost
  estimate. Heuristics: rounding, a QP dive and feasibility jump on the objective linearized
  at the root solution. The barrier remains the fallback when the active-set method fails.
  Compared with barrier relaxations at every node: miqprand_large 4.9 s -> 0.3 s,
  miqprand_hard (new) 43 s -> 2.7 s, miqpport_hard (new) 13.4 s -> 0.8 s.
- C API: `popso_set_q`, `popso_add_q_entries`, `popso_num_q_nonzeros`, `popso_get_q`. Python:
  `popso.QuadExpr` (products and squares of variables and linear expressions) in
  `set_objective`, `Model.get_q()`.
- Benchmarks: new families `portfolio`, `lsq` (constrained least squares) and the MIQPs
  `miqpport` (cardinality-constrained portfolio) and `miqprand`; the solution checker
  recomputes quadratic objectives; `benchmarks/tools/qp_fuzz.py` compares random QPs and MIQPs
  with HiGHS and SCIP; `benchmarks/tools/maros_meszaros.py` fetches and converts the
  Maros–Meszaros QP test set. The 16 generated QP/MIQP instances are solved in a shifted
  geometric mean of 0.25 s (2 threads; references: HiGHS 1.5 s on 14, SCIP 5.6 s on 12).
- Maros–Meszaros (qpbenchmark's collection, 2 threads): the 98 problems up to 60 KB are all
  solved (97 convex, 1.8 s in total; VALUES is reported as nonconvex); of the 40 larger ones
  (60 s limit) 32 are solved (BOYD1/2, CONT-300, EXDATA among them); LISWET1, 7-10 and 12 end
  with a numerical error and CVXQP1_L / CVXQP3_L at the time limit.

**Symmetry**
- Symmetry-breaking reductions from a stabilizer chain (Schreier-Sims table / SST): after
  presolve, leaders whose orbit contains columns in conflict with them (a row forbids both at
  1) give x_leader >= x_follower rows and fix the conflicting followers to 0, then the search
  continues in (a subgroup of) the leader's stabilizer. Used when the fixings remove at least
  10% of the columns; the problem is presolved again and the remaining symmetry is handled as
  before (packing orbitopes, orbital branching). cod105: 1024 -> 341 columns, solved at the
  root (12 s -> 1-3 s).
- CLI `--symmetry on|off` (default on) switches all symmetry handling.

**Fixes**
- Barrier: inequality rows with fixed columns had their slack bounds shifted twice (the
  crossover repaired the point, so LP results were right, but the interior solution was not).

**Cutting planes**
- c-MIR path aggregation now eliminates a continuous column through the tight row with the largest
  LP dual (the rows that prove the LP bound) instead of the first tight row found (not for columns
  with a simple variable upper bound x <= u y, where the old order suits fixed-charge paths). On the
  63-instance MIPLIB subset the root node closes 65.2% of the root gap on average (63.8% before;
  SCIP 64.6% with node limit 1); neos17 goes from 18% to 94% of its root gap closed and solves in
  about 1 s instead of 4-6 s. End to end the subset is unchanged (60 s / 120 s A/Bs within noise).

## v0.3.1 — 2026-11

Benchmark update; the solver is the same as in v0.3.0.

- BENCHMARKS.md: new section on 43 hard MIPLIB 2017 instances (those of the 106-instance
  development set that no solver finishes within 120 s), run with a 10-minute limit. popso solves
  9, HiGHS 12, SCIP 6; set file `benchmarks/sets/miplib_hard.txt`.
- HiGHS and SCIP reference results for these instances at 600 s.
- `plot_profile.py`: `--tmin` and `--ymax` to zoom the time and count axes.

## v0.3.0 — 2026-11

New capabilities: an interior-point (barrier) LP solver and the logical modelling features real
models use — semi-continuous variables, SOS1/SOS2 and indicator constraints.

**Benchmarks**
- 63-instance MIPLIB 2017 subset (120 s, 4 threads): 60 solved, shifted geometric mean 11.4 s
  (unchanged from v0.2.1: this release adds features, not MILP speed).
- Generated models, now 130 including three families with logical constraints (1 thread):
  129 solved; SCIP solves 122 of 130, HiGHS 117 of the 122 it can read.
- LP relaxations of the MIPLIB instances: the new default `--lp-method auto` is about 27% faster
  than the dual simplex alone. BENCHMARKS.md has a new plot of popso's progress across versions.

**Barrier (interior-point) LP solver**
- New LP algorithm, written from scratch like the rest: Mehrotra's predictor-corrector
  primal-dual method with Gondzio's centrality correctors, on the presolved and scaled problem
  (geometric scaling; free variables, ranges and two-sided bounds handled directly).
- Linear algebra: an approximate minimum degree ordering (AMD) and a supernodal sparse Cholesky
  factorization, multithreaded for large dense blocks. Each solve picks the cheaper of the normal
  equations and the regularized augmented system (factored as a quasi-definite LDL^T, which suits
  LPs with far more rows than columns and is the form a future QP solver needs). Dense columns in
  the normal equations go through a Sherman–Morrison–Woodbury update; every solve is refined
  iteratively.
- Crossover to an optimal basis: dual and primal push phases, then a short simplex cleanup, so
  the answer is an exact vertex as with the simplex.
- `--lp-method auto|simplex|barrier` (C API parameter `lp_method`: 0, 1, 2; the Python bindings
  also take the names). The method is used for a pure LP and for the root LP of a MIP. `auto`
  (the default) runs the dual simplex and the barrier at the same time on large pure LPs and
  takes the first to finish (not in `--deterministic` mode, where it stays with the simplex);
  MIP root LPs stay with the simplex unless `barrier` is chosen.
- On LP relaxations of 69 MIPLIB instances (4 threads), `auto` takes 27% less time than the
  dual simplex alone in shifted geometric mean and 44% less in total; the barrier alone is about
  as fast as HiGHS' interior-point solver (shifted geometric mean 0.40 s vs 0.43 s, with
  crossover). Examples: chromaticindex512-7 20 s -> 1.6 s, dano3_3 6.2 s -> 2.7 s, qap10
  3.4 s -> 1.1 s, rmatr200-p5 5.9 s -> 2.9 s, brazil3 2.4 s -> 0.8 s.

**Model features**
- Semi-continuous and semi-integer variables (x = 0 or lb <= x <= ub), SOS1 and SOS2
  constraints, and indicator constraints (if a binary variable takes a given value, a linear
  constraint must hold).
- They are read from MPS files (`SC`/`SI` bounds, `SOS` section, `INDICATORS` section) and LP
  files (`Semi-Continuous` section, `SOS` section with `S1::`/`S2::` entries, and indicator
  constraints written `name: b = 1 -> x + y <= 5`). The MPS writer writes them, so models
  round-trip.
- C API: `POPSO_SEMICONTINUOUS` and `POPSO_SEMIINTEGER` variable types, `popso_add_sos`,
  `popso_add_indicator`, `popso_num_sos`, `popso_num_indicators`, and `popso_read_lp` for LP
  files. Python: `SEMICONTINUOUS`/`SEMIINTEGER` variable types, `Model.add_sos`,
  `Model.add_indicator(b, 1, x + y <= 5)` and `popso.read_lp`.
- How they are solved: the solver adds a binary variable for each semi-continuous variable, SOS
  member or SOS2 segment and links it with big-M rows; an indicator constraint becomes a big-M
  row. The big-M values come from the variable bounds, tightened by bound propagation over the
  model rows. Presolve, cuts, heuristics and branching then work as for any MILP, and the
  final solution is checked against the original constraints. When no finite bound exists (or
  it would exceed 1e6), the constraint is not relaxed by a big-M: the solver branches on its
  binary variable and enforces the constraint in that subtree. Presolve is skipped for such
  models.
- Testing: unit tests compare small models with optima found by enumeration, using both the
  big-M form and branching only, and `benchmarks/tools/logical_fuzz.py` compares random small
  models with SCIP. Three new generated benchmark families, `semilot` (minimum
  lot sizes), `pwlnet` (transportation with piecewise-linear concave costs, plant opening and
  single sourcing) and `indsched` (single-machine scheduling with indicator disjunctions and no
  big-M), have reference answers from SCIP (and from HiGHS for `semilot`). The benchmark
  checker verifies the new constraints.

## v0.2.1 — 2026-11

A patch release: faster root cuts on pure integer models and a robustness fix for deterministic
parallel mode under tight memory limits.

**Benchmarks**
- 63-instance MIPLIB 2017 subset (120 s, 4 threads): 60 solved, shifted geometric mean 11.3 s
  (v0.2.0: 12.4 s). See [BENCHMARKS.md](BENCHMARKS.md).

**Solver**
- Deterministic parallel mode (`--deterministic`) now drains the search tree reliably near the
  memory cap (`--mem-limit`). Before, a parallel deterministic run could stop with `MEMORY_LIMIT`
  where a 1-thread run under the same cap proved optimality. The work between two
  synchronization points is now sized so that the tree stays below the cap. Reruns stay
  identical, and runs that stay well below the cap are unchanged.

**Root cuts**
- c-MIR path aggregation continues through general integer columns on pure integer models (root cut loop,
  work-capped, not on dual degenerate root LPs): mzzv11 / mzzv42z are now solved at the root (4 threads:
  51 s -> 17 s and 15 s -> 3.3 s).
- Clique cuts are extended to maximal cliques by zero-valued literals.
- `benchmarks/tools/rootgap.py`: root gap closed per instance vs SCIP / HiGHS (node limit 1) and the known
  optimum; the solver prints a machine-readable `ROOTSTAT` line after each root.

## v0.2.0 — 2026-11

popso becomes usable as a library and on all three major platforms: a C API and Python bindings,
starting solutions, an LP-format reader, CI on Linux, macOS and Windows, a stronger presolve and a
more robust LP factorization.

**Benchmarks**
- 63-instance MIPLIB 2017 subset (120 s, 4 threads): 60 solved (v0.1.0: 59), shifted geometric
  mean 12.4 s (v0.1.0: 13.2 s); HiGHS solves 53, SCIP 50. 118 generated models (1 thread): 116
  solved. See [BENCHMARKS.md](BENCHMARKS.md).

**Solver**
- Stronger presolve:
  - clique extension: set-packing rows grow with conflicting binaries, dominated packing rows are
    dropped, and binaries that conflict with a whole equality row are fixed;
  - dual bound strengthening and GCD-based tightening of integer rows;
  - a set-cover rule for dominated columns;
  - faster probing that also tightens continuous bounds.

  On the MIPLIB subset the presolved models now close more of the root gap on average than
  SCIP's or HiGHS's presolve, and presolve takes about 30% less time than in v0.1.0.
- LP factorization: an update that would make the factorization inaccurate (a very large row
  multiplier in the Forrest–Tomlin update) now triggers a refactorization instead of carrying the
  error forward.

**Library API**
- C API in `include/popso.h`, built as the shared library `popso` (`libpopso.so`,
  `libpopso.dylib`, `popso.dll`). It can build models (variables, sparse rows, objective) or read
  MPS files, and it sets parameters, solves, and queries status, objective, bound, gap, node count,
  time and the solution. It can also write MPS and solution files. Calls return error codes, and
  `popso_last_error` gives a message. Separate models can be solved from different threads.
- Python bindings (`python/popso`), using `ctypes` with no compiled extension and no
  dependencies. They support linear expressions (`m.add_constr(x + 2 * y <= 5)`), a result object
  that can be indexed by variable or name, parameters, and starting solutions. The package can be
  installed with pip.

**Starting solutions**
- `--start <file>` reads a starting solution in the `name value` format that `--sol` writes. A
  complete feasible start becomes the first incumbent. A partial or infeasible start keeps its
  integer values and is completed by a short sub-MIP; if that fails it is ignored, and the log
  says so. Starts are also available in the C API and in Python.

**Input / output**
- CPLEX LP-format reader (`.lp` files, as written by CPLEX, HiGHS, SCIP and modelling
  tools): objective constant, ranged rows (`lo <= expr <= hi`), all bound forms
  (`free`, `-inf`/`+inf`, fixed, double-sided), General/Binary sections, comments and
  case-insensitive keywords. Quadratic terms, indicator constraints, semi-continuous
  variables and SOS are rejected with a clear "not supported yet" error and line number.
- The command line picks the reader by file extension (`.lp` vs `.mps`, case-insensitive);
  `--format mps|lp` overrides it.
- Testing: a pure-Python LP writer in the benchmark modeller and
  `benchmarks/tools/lp_crosscheck.py`, which checks that every benchmark instance read as
  LP gives exactly the same model as the MPS file, and solves LP files written by HiGHS
  and SCIP.

**Platforms**
- Continuous integration on GitHub Actions: every push and pull request is built and tested on
  Linux (gcc and clang), macOS (clang) and Windows (MSVC). Each build runs the unit tests, the
  C API and Python binding tests, and solves the tiny benchmark instances with independently
  verified answers. A separate Linux build runs the unit and C API tests under
  AddressSanitizer and UndefinedBehaviorSanitizer.
- Builds without warnings with gcc 15 and clang. The new CMake option `POPSO_WERROR` turns
  warnings into errors.
- Builds with gcc on macOS (the macOS system headers used for memory measurement did not
  compile with gcc).
- Windows/MSVC: `windows.h` no longer defines `min`/`max` macros, there are no warnings about
  standard C functions such as `getenv`, and the MSVC warning set matches the gcc/clang one.
- Benchmark runner: `benchmarks/run.py` needs only the Python standard library unless it runs
  the reference solvers. It finds `popso.exe` and `build/Release/` on Windows, works there
  without `ps`, and has a `--strict` flag that fails when an instance is not solved.

## v0.1.0 — 2026-11

First public release: a complete, multithreaded LP/MILP solver written from scratch in
C++20, using only the standard library.

**Solves**
- Linear programs (LP) and mixed-integer linear programs (MILP) with continuous, integer
  and binary variables, ranged rows and free variables.
- It proves optimality, infeasibility or unboundedness.
- On a 63-instance MIPLIB 2017 subset it solves 59 within 120 s on 4 threads, against 53
  for HiGHS and 50 for SCIP. On 118 generated models it solves 116 on 1 thread, against
  113 for HiGHS and 111 for SCIP. See [BENCHMARKS.md](BENCHMARKS.md).

**Input / output**
- MPS reader (free and fixed format, including RANGES, OBJSENSE and integer markers).
- Solution file output (`name value` per line), plus a summary of status, objective, bound,
  gap, nodes and time.
- Command-line interface with time, gap, node and memory limits, thread count, seed,
  deterministic mode and racing.

**LP engine**
- Bounded dual simplex.
- LU factorization with Forrest–Tomlin updates and hyper-sparse solves.
- Dual steepest-edge pricing, cost perturbation, and a primal simplex for cleanup.
- Warm starts from stored bases throughout the tree.

**Presolve**
- Removes singleton rows and columns, parallel rows and columns, and dominated columns.
- Tightens bounds and coefficients.
- Aggregates and substitutes variables.
- Probing with implications, and clique merging.
- Full postsolve; restarts after root fixings.

**Cutting planes**
- Gomory mixed-integer cuts and c-MIR cuts on aggregated rows (with variable-bound
  substitution).
- Lifted knapsack covers, also on aggregated rows.
- Clique, implied-bound and zero-half cuts.
- Tableau cuts at tree nodes.
- A dynamic cut pool with aging.

**Primal heuristics**
- Rounding, diving, local search.
- Feasibility pump and feasibility jump.
- Large-neighbourhood sub-MIPs (RINS, RENS, crossover) with a solution pool.

**Branch-and-bound**
- Branching: reliability branching on pseudocosts with strong branching.
- Pruning:
  - node-level domain propagation;
  - conflict analysis;
  - reduced-cost fixing.
- Symmetry: orbital fixing and packing orbitopes.
- Restarts at the root and in the tree.
- Node selection mixes best bound and best estimate.
- Recovery from numerical failures in node LPs.

**Parallelism and resources**
- Multithreaded tree search, opportunistic by default.
- A deterministic mode that gives the same result on every run with the same thread count.
- Search-tree memory stays within a configurable cap, using compact node storage; near the
  cap the solver finishes off open nodes rather than creating new ones.

**Known limitations**
- No quadratic or nonlinear models.
- No SOS or indicator constraints, and no semi-continuous variables.
- No LP-format reader, no library API, and no barrier/interior-point LP method.
