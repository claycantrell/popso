#pragma once
// Improvement heuristics for a feasible point: 1-opt shifting and a local
// branching sub-MIP (Fischetti-Lodi) around an incumbent.

#include <vector>

#include "milp/heuristics/submip.hpp"
#include "milp/search/mip_problem.hpp"

namespace milp {

// Shifts single variables in their objective-improving direction as far as
// the rows and bounds allow (integers by whole steps). `x` must be feasible;
// it stays feasible. Returns true if the objective improved.
bool one_opt(const MipProblem& p, std::vector<double>& x);

// Solves p with the extra row  sum_{binaries} |x_j - inc_j| <= k.
// Returns true with a solution whose objective (incl. offset) is below `cutoff`.
bool local_branching(const MipProblem& p, const std::vector<double>& incumbent, int k, double cutoff,
                     const SubMipLimits& limits, std::vector<double>& sol);

}  // namespace milp
