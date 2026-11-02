#pragma once
// Probing on binary variables: each binary x is tentatively fixed to 0 and to
// 1 and the rows are propagated (activity-based).
//   - a setting that is infeasible fixes x to the other value;
//   - bounds implied by both settings hold globally;
//   - a column y fixed to v0 by x = 0 and to v1 by x = 1 is y = v0 + (v1 - v0) x
//     (for binaries: y = x or y = 1 - x).
// Every derivation is sound with conservative tolerances, so a work limit only
// makes the result weaker.

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
  long long probes = 0, work = 0;
};

// work_limit: roughly the number of matrix entries visited.
ProbingResult probe_binaries(const MipProblem& p, long long work_limit);

}  // namespace milp
