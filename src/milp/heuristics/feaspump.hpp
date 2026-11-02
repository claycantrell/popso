#pragma once
// Objective feasibility pump (Fischetti-Glover-Lodi; Achterberg-Berthold):
// alternates between rounding the LP point and projecting the rounding back
// onto the LP polyhedron (L1 distance over integer variables), with random
// perturbation to escape cycles.

#include <vector>

#include "milp/lp/simplex.hpp"
#include "milp/search/mip_problem.hpp"

namespace milp {

struct FeasPumpLimits {
  int max_rounds = 200;
  double time_limit = 10.0;          // seconds
  long long iteration_budget = 100000;  // simplex iterations
  unsigned seed = 0;
  int stall_rounds = 20;  // stop after this many rounds without a closer projection (0: no limit)
};

// `root` holds the solved root LP (its basis is used for warm starts).
// Returns true with a feasible point in `sol`.
// On failure, `closest` (if given) receives the LP point that was closest to
// its rounding, a good centre for a follow-up sub-MIP (FP "stage 3").
bool feasibility_pump(const MipProblem& p, const DualSimplex& root, const FeasPumpLimits& limits,
                      std::vector<double>& sol, long long* iterations, std::vector<double>* closest = nullptr);

}  // namespace milp
