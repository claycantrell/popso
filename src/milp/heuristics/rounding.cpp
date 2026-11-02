#include "milp/heuristics/rounding.hpp"

#include <cmath>

namespace milp {

SimpleRounding::SimpleRounding(const MipProblem& p) : p_(p) {
  const Index n = p.num_cols();
  down_locks_.assign(n, 0);
  up_locks_.assign(n, 0);
  for (Index j = 0; j < n; ++j) {
    for (Index k = p.A.begin(j); k < p.A.end(j); ++k) {
      const Index i = p.A.index[k];
      const bool has_lb = is_finite(p.row_lb[i]), has_ub = is_finite(p.row_ub[i]);
      if (p.A.value[k] > 0) {
        down_locks_[j] += has_lb;
        up_locks_[j] += has_ub;
      } else {
        down_locks_[j] += has_ub;
        up_locks_[j] += has_lb;
      }
    }
  }
}

bool SimpleRounding::run(const double* x, std::vector<double>& out) const {
  const Index n = p_.num_cols();
  out.assign(x, x + n);
  for (Index j = 0; j < n; ++j) {
    if (!p_.is_int[j]) continue;
    const double v = x[j];
    const double r = std::round(v);
    if (std::abs(v - r) <= 1e-6) {
      out[j] = r;
      continue;
    }
    if (down_locks_[j] == 0) out[j] = std::floor(v);
    else if (up_locks_[j] == 0) out[j] = std::ceil(v);
    else return false;
  }
  return p_.max_violation(out) <= 1e-6;
}

}  // namespace milp
