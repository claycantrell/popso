# Changelog

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
