# LP engine stress tests

- `lu_random.cpp` — random sparse/singular bases: factorize, eta updates, ftran/btran residuals (dense and
  HVector solves), plus large LP-like bases (m up to 1500). `./lu_random [updates]`; `FRESH=1` refactors after
  the updates (separates factorization error from eta-update error); `THR=<u>` sets the pivot threshold.
- `lp_random.cpp` — random LPs (all bound types, ranged rows): KKT check of optimal answers, and
  dumps every LP + our status to `lps.txt`. `./lp_random <trials> [size]`: size 10 gives LPs up to
  400 x 600 that exercise the Markowitz kernel and dense LU block (HiGHS may report 'unknown' on some).
- `dive_bench.cpp` — LP engine benchmark and warm-start checker on a MIP model: simulates the node LP sequence of
  a depth-first search (bound changes, warm solves, backtracking via `restore_state`/`set_basis`, iteration-limited
  probes) and reports us/solve and us/iteration. `./dive_bench model.mps [solves] [check_every]`; `check_every=k`
  compares every k-th warm solve with a cold solve (status + objective).
- `MILP_LU_STATS=1` (any program) — per-factorization stats on stderr and a total at exit.
- `MILP_LU_HYPER=1` / `=2` (any program) — LU solves always use the dense loops / always the hyper-sparse
  path (default: hyper-sparse when the right-hand side and the predicted result are below 10% density).
- `xcheck.py lps.txt` — re-solves every dumped LP with HiGHS and reports status/objective mismatches.

```
c++ -std=c++20 -O2 -Isrc tools/lpcheck/lp_random.cpp src/milp/lp/*.cpp src/milp/core/*.cpp -o lp_random
./lp_random 3000 && .venv/bin/python tools/lpcheck/xcheck.py lps.txt
c++ -std=c++20 -O2 -Isrc tools/lpcheck/dive_bench.cpp src/milp/lp/*.cpp src/milp/core/*.cpp src/milp/io/mps_reader.cpp -o dive_bench
```

# Search debugging switches (environment variables)

- `MILP_DEBUG_SOL=<file>` — file with a known optimal point (one value per line, column order).
  Logs every node whose bounds contain that point when it gets closed, and gap-closure events.
  HiGHS can produce such a file: `h.getSolution().col_value`. Also reports every learned conflict
  constraint that the point violates although it is better than the cutoff the conflict relies on.
- `MILP_CHECK_LP=1` — re-solves every node LP from scratch and reports warm/fresh mismatches.
- `MILP_NO_SB=1`, `MILP_NO_DIVE=1` — disable strong branching / diving heuristics.
- `MILP_NO_CONFLICT=1` — disable conflict analysis (dual proofs of infeasible / cut-off nodes).
