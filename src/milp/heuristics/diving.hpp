#pragma once
// LP diving heuristics: repeatedly round one fractional integer variable
// (by tightening its bound) and resolve the LP until the solution is
// integral, the LP becomes infeasible, or the budget runs out.

#include <vector>

#include "milp/heuristics/rounding.hpp"
#include "milp/lp/simplex.hpp"
#include "milp/search/mip_problem.hpp"
#include "milp/search/propagation.hpp"

namespace milp {

enum class DiveRule { Fractional, Coefficient, Guided, VectorLength };
const char* to_string(DiveRule r);

struct DiveLimits {
  long long iteration_budget = 5000;  // total simplex iterations for the dive
  double time_limit = kInf;           // seconds
  double cutoff = kInf;               // abandon once the LP objective reaches this
};

// `lp` must hold an optimal LP solution of the starting node. On return the
// LP bounds are restored, but its basis is left wherever the dive ended
// (callers re-install the node basis). `incumbent` may be null (Guided then
// behaves like Fractional). With `prop`, each rounding is propagated and the
// implied integer bounds are applied to the LP too. Returns true with a
// feasible point in `sol`.
bool dive(const MipProblem& p, DualSimplex& lp, DiveRule rule, const SimpleRounding& rounding,
          const std::vector<double>* incumbent, const DiveLimits& limits, std::vector<double>& sol,
          long long* iterations, DomainPropagator* prop = nullptr);

}  // namespace milp
