#include "milp/presolve/implied_int.hpp"

#include <cmath>
#include <vector>

namespace milp {

namespace {

bool integral(double v) { return std::abs(v - std::round(v)) <= 1e-9 * std::max(1.0, std::abs(v)); }

}  // namespace

std::vector<char> detect_implied_integers(const MipProblem& p, bool dual, double bound_tol, ImpliedIntStats* stats) {
  const Index n = p.num_cols(), m = p.num_rows();
  std::vector<char> impl(static_cast<size_t>(n), 0);
  if (n == 0) return impl;
  dual = dual && !p.quadratic() && p.implications.empty();
  const SparseMatrix R = p.A.transpose();  // CSR: row i's entries in R.begin(i) .. R.end(i)
  // intlike: integer, integrally fixed, or implied integer (found so far).
  std::vector<char> intlike(static_cast<size_t>(n), 0);
  for (Index j = 0; j < n; ++j) {
    if (p.is_int[j]) {
      intlike[j] = 1;
    } else if (p.is_implied_int(j)) {
      intlike[j] = impl[j] = 1;
    } else if (p.col_lb[j] == p.col_ub[j] && integral(p.col_lb[j])) {
      intlike[j] = 1;
    }
  }
  std::vector<Index> nonint(static_cast<size_t>(m), 0);
  for (Index i = 0; i < m; ++i) {
    for (Index k = R.begin(i); k < R.end(i); ++k) nonint[i] += !intlike[R.index[k]];
  }
  // For a row with exactly one non-intlike column j: whether the row makes x_j integral (equation with
  // integral ratios and b / a_j) and whether it allows the dual rule (integral ratios and sides / a_j).
  std::vector<char> row_eq(static_cast<size_t>(m), 0), row_dual(static_cast<size_t>(m), 0);
  std::vector<char> row_done(static_cast<size_t>(m), 0);
  auto eval_row = [&](Index i) {
    row_done[i] = 1;
    Index j = -1;
    double aj = 0.0;
    for (Index k = R.begin(i); k < R.end(i); ++k) {
      if (!intlike[R.index[k]]) j = R.index[k], aj = R.value[k];
    }
    if (j < 0 || aj == 0.0) return;
    bool ratios = true;
    for (Index k = R.begin(i); k < R.end(i) && ratios; ++k) {
      if (R.index[k] != j) ratios = integral(R.value[k] / aj);
    }
    if (!ratios) return;
    const double lo = p.row_lb[i], hi = p.row_ub[i];
    const bool lo_ok = !std::isfinite(lo) || integral(lo / aj), hi_ok = !std::isfinite(hi) || integral(hi / aj);
    row_eq[i] = lo == hi && std::isfinite(lo) && lo_ok;
    row_dual[i] = lo_ok && hi_ok;
  };
  auto bound_ok = [&](double b) { return !std::isfinite(b) || std::abs(b - std::round(b)) <= bound_tol; };
  // Worklist of columns to (re)check: initially every non-intlike column; a column is queued again when
  // one of its rows drops to a single non-intlike column.
  std::vector<Index> work;
  std::vector<char> queued(static_cast<size_t>(n), 0);
  for (Index j = 0; j < n; ++j) {
    if (!intlike[j] && p.A.end(j) > p.A.begin(j)) work.push_back(j), queued[j] = 1;
  }
  Index n_primal = 0, n_dual = 0;
  while (!work.empty()) {
    const Index j = work.back();
    work.pop_back();
    queued[j] = 0;
    if (intlike[j]) continue;
    bool primal = false, dual_ok = dual && bound_ok(p.col_lb[j]) && bound_ok(p.col_ub[j]);
    for (Index k = p.A.begin(j); k < p.A.end(j); ++k) {
      const Index i = p.A.index[k];
      if (p.A.value[k] == 0.0) continue;
      if (!std::isfinite(p.row_lb[i]) && !std::isfinite(p.row_ub[i])) continue;  // free row
      if (nonint[i] != 1) {
        dual_ok = false;
        continue;
      }
      if (!row_done[i]) eval_row(i);
      primal = primal || row_eq[i];
      dual_ok = dual_ok && row_dual[i];
    }
    if (!primal && !dual_ok) continue;
    intlike[j] = impl[j] = 1;
    ++(primal ? n_primal : n_dual);
    for (Index k = p.A.begin(j); k < p.A.end(j); ++k) {
      const Index i = p.A.index[k];
      if (--nonint[i] != 1) continue;
      for (Index q = R.begin(i); q < R.end(i); ++q) {
        const Index c = R.index[q];
        if (!intlike[c] && !queued[c]) work.push_back(c), queued[c] = 1;
      }
    }
  }
  if (stats) stats->primal += n_primal, stats->dual += n_dual;
  return impl;
}

Index mark_implied_integers(MipProblem& p, bool dual, ImpliedIntStats* stats) {
  std::vector<char> impl = detect_implied_integers(p, dual, 1e-6, stats);
  Index cnt = 0;
  for (Index j = 0; j < p.num_cols(); ++j) {
    if (!impl[j]) continue;
    double l = p.col_lb[j], u = p.col_ub[j];
    if (std::isfinite(l)) l = std::ceil(l - 1e-6);
    if (std::isfinite(u)) u = std::floor(u + 1e-6);
    if (l > u) {
      impl[j] = 0;  // no integral value: the problem is infeasible; leave it to the search
      continue;
    }
    p.col_lb[j] = l;
    p.col_ub[j] = u;
    ++cnt;
  }
  p.implied_int = cnt > 0 ? std::move(impl) : std::vector<char>{};
  return cnt;
}

}  // namespace milp
