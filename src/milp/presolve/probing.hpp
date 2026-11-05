#pragma once
// Probing on binary variables: each binary x is tentatively fixed to 0 and to
// 1 and the rows are propagated (activity-based).
//   - a setting that is infeasible fixes x to the other value;
//   - bounds implied by both settings hold globally;
//   - a column y fixed to v0 by x = 0 and to v1 by x = 1 is y = v0 + (v1 - v0) x
//     (for binaries: y = x or y = 1 - x);
//   - lifting: a one-sided row whose maximal activity (<= form, side b) under the
//     bounds implied by x = v is b - delta < b gets x into it: v = 1: a_x += delta
//     (x = 1 then leaves the row the slack it never used), v = 0: a_x -= delta and
//     b -= delta. The lifted row cuts off no integer-feasible point and implies the
//     original one for x in [0, 1];
//   - implications between binaries, x = v implies y = w, kept as conflicting literals
//     (x = v, y = 1 - w) and propagated together with the rows by later probes (so
//     also backwards: y = 1 - w implies x = 1 - v, which the rows alone may not
//     show); returned for later passes and, x + y <= 1, for the clique extension of
//     set-packing rows (presolve/clique_merge.hpp).
// Every derivation is sound with conservative tolerances, so a work limit only
// makes the result weaker.

#include <utility>
#include <vector>

#include "milp/core/types.hpp"
#include "milp/search/mip_problem.hpp"

namespace milp {

struct ProbingResult {
  bool infeasible = false;
  std::vector<double> lb, ub;  // tightened global bounds
  struct Equiv {
    Index y, x;
    double v0, v1;  // y = v0 + (v1 - v0) x
  };
  std::vector<Equiv> equiv;
  // Liftings, at most one per row (two liftings of the same row are not jointly valid: each assumes the
  // original row). upper: the row_ub side (<= form); else the row_lb side, in <= form -row <= -row_lb.
  // In <= form: one (v = 1): a_col += delta; else a_col -= delta and the side -= delta.
  struct Lift {
    Index row, col;
    double delta;
    bool upper, one;
    bool fill;     // col is not in the row (a new entry)
    double score;  // delta / max |a| of the row
  };
  std::vector<Lift> lifts;
  // New conflicting literal pairs a < b (literal 2 j + v: x_j = v; a and b are not both true).
  std::vector<std::pair<Index, Index>> conflicts;
  long long probes = 0, work = 0;
};

struct ProbingOptions {
  bool lift = false;         // find liftings (the rows of p are not changed; the caller applies them)
  size_t max_conflicts = 0;  // keep and propagate up to this many conflicting literal pairs (incl. seeds)
  const std::vector<std::pair<Index, Index>>* seeds = nullptr;  // known conflicts (binaries of p)
};

// work_limit: roughly the number of matrix entries visited.
ProbingResult probe_binaries(const MipProblem& p, long long work_limit, const ProbingOptions& opt = {});

}  // namespace milp
