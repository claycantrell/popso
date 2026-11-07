#include "milp/heuristics/rounding.hpp"

#include <algorithm>
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

bool ZiRounding::run(const double* x, std::vector<double>& out, long long* work) const {
  const Index n = p_.num_cols(), m = p_.num_rows();
  const SparseMatrix& A = p_.A;
  constexpr double kFrac = 1e-6;
  std::vector<Index> frac;
  for (Index j = 0; j < n; ++j) {
    if (p_.is_int[j] && std::abs(x[j] - std::round(x[j])) > kFrac) frac.push_back(j);
  }
  if (frac.empty()) return false;
  out.assign(x, x + n);
  std::vector<double> act(m, 0.0);
  for (Index j = 0; j < n; ++j) {
    if (out[j] == 0.0) continue;
    for (Index k = A.begin(j); k < A.end(j); ++k) act[A.index[k]] += A.value[k] * out[j];
  }
  long long touched = static_cast<long long>(A.nnz());
  // Room to move column j up / down without leaving its bounds or any row's sides (a side the point already
  // violates within the LP tolerance gives no room).
  auto room = [&](Index j, double* up, double* down) {
    *up = p_.col_ub[j] - out[j];
    *down = out[j] - p_.col_lb[j];
    for (Index k = A.begin(j); k < A.end(j); ++k) {
      const Index i = A.index[k];
      const double a = A.value[k];
      const double su = std::max(0.0, p_.row_ub[i] - act[i]), sd = std::max(0.0, act[i] - p_.row_lb[i]);
      if (a > 0.0) {
        *up = std::min(*up, su / a);
        *down = std::min(*down, sd / a);
      } else {
        *up = std::min(*up, sd / -a);
        *down = std::min(*down, su / -a);
      }
    }
    touched += A.end(j) - A.begin(j);
  };
  for (int pass = 0; pass < 3 && !frac.empty(); ++pass) {
    size_t keep = 0;
    for (Index j : frac) {
      const double v = out[j], lo = std::floor(v), hi = std::ceil(v);
      double up, down;
      room(j, &up, &down);
      const bool can_down = v - lo <= down + 1e-9, can_up = hi - v <= up + 1e-9;
      double to;
      if (can_down && can_up) {
        // Both: the cheaper direction, else the nearer integer.
        const double c = p_.cost[j];
        to = c > 0.0 ? lo : c < 0.0 ? hi : (v - lo <= hi - v ? lo : hi);
      } else if (can_down) {
        to = lo;
      } else if (can_up) {
        to = hi;
      } else {
        frac[keep++] = j;
        continue;
      }
      const double d = to - v;
      out[j] = to;
      for (Index k = A.begin(j); k < A.end(j); ++k) act[A.index[k]] += A.value[k] * d;
    }
    if (keep == frac.size()) break;  // no shift in this pass
    frac.resize(keep);
  }
  if (work) *work += touched;
  if (!frac.empty()) return false;
  return p_.max_violation(out) <= 1e-6;
}

}  // namespace milp
