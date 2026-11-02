#pragma once
// Lock-based simple rounding: rounds each fractional integer variable in a
// direction in which no row can become violated. Succeeds only if every
// fractional variable has such a direction (then feasibility is guaranteed).

#include <vector>

#include "milp/search/mip_problem.hpp"

namespace milp {

class SimpleRounding {
 public:
  explicit SimpleRounding(const MipProblem& p);
  // x: LP solution. On success writes a feasible point to `out`.
  bool run(const double* x, std::vector<double>& out) const;

  Index down_locks(Index j) const { return down_locks_[j]; }
  Index up_locks(Index j) const { return up_locks_[j]; }

 private:
  const MipProblem& p_;
  std::vector<Index> down_locks_, up_locks_;
};

}  // namespace milp
