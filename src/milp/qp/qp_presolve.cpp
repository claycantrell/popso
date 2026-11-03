#include "milp/qp/qp_presolve.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>

#include "milp/qp/qp.hpp"

namespace milp {

namespace {

constexpr double kFeasTol = 1e-6;      // infeasibility proofs (relative to 1 + |bound|)
constexpr double kRedundantTol = 1e-10;  // a row side is dropped if the activity bound reaches it this closely
constexpr int kMaxRounds = 25;

double tol_of(double v) { return kFeasTol * (1.0 + std::abs(v)); }

}  // namespace

QpPresolveStatus QpPresolver::run(const MipProblem& p, const Logger& log) {
  orig_ = p;
  const Index n = p.num_cols(), m = p.num_rows();
  stats_ = QpPresolveStats();
  stats_.rows_before = m;
  stats_.cols_before = n;
  std::vector<double> lb = p.col_lb, ub = p.col_ub, c = p.cost, rlo = p.row_lb, rhi = p.row_ub;
  double offset = p.offset;
  const bool has_int = std::any_of(p.is_int.begin(), p.is_int.end(), [](char v) { return v != 0; });

  // Column- and row-wise copies of A; Hessian: diagonal and symmetric adjacency.
  std::vector<std::vector<std::pair<Index, double>>> col_a(static_cast<size_t>(n)), row_a(static_cast<size_t>(m));
  for (Index j = 0; j < n; ++j) {
    for (Index k = p.A.begin(j); k < p.A.end(j); ++k) {
      if (p.A.value[k] == 0.0) continue;
      col_a[j].emplace_back(p.A.index[k], p.A.value[k]);
      row_a[p.A.index[k]].emplace_back(j, p.A.value[k]);
    }
  }
  std::vector<double> qdiag(static_cast<size_t>(n), 0.0);
  std::vector<std::vector<std::pair<Index, double>>> qadj(static_cast<size_t>(n));
  for (Index j = 0; j < p.Q.num_major; ++j) {
    for (Index k = p.Q.begin(j); k < p.Q.end(j); ++k) {
      const Index i = p.Q.index[k];
      const double v = p.Q.value[k];
      if (v == 0.0) continue;
      if (i == j) {
        qdiag[j] += v;
      } else {
        qadj[j].emplace_back(i, v);
        qadj[i].emplace_back(j, v);
      }
    }
  }
  std::vector<Index> row_len(static_cast<size_t>(m), 0), col_len(static_cast<size_t>(n), 0), q_len(static_cast<size_t>(n), 0);
  for (Index j = 0; j < n; ++j) {
    col_len[j] = static_cast<Index>(col_a[j].size());
    q_len[j] = static_cast<Index>(qadj[j].size());
  }
  for (Index i = 0; i < m; ++i) row_len[i] = static_cast<Index>(row_a[i].size());
  col_removed_.assign(static_cast<size_t>(n), 0);
  std::vector<char> row_removed(static_cast<size_t>(m), 0);
  fixed_val_.assign(static_cast<size_t>(n), 0.0);
  lb_row_.assign(static_cast<size_t>(n), -1);
  ub_row_.assign(static_cast<size_t>(n), -1);
  singletons_.clear();
  bool infeasible = false;

  auto fix_col = [&](Index j, double v) {
    offset += c[j] * v + 0.5 * qdiag[j] * v * v;
    for (const auto& [k, q] : qadj[j]) {
      if (col_removed_[k]) continue;
      c[k] += q * v;
      --q_len[k];
    }
    for (const auto& [i, a] : col_a[j]) {
      if (row_removed[i]) continue;
      if (is_finite(rlo[i])) rlo[i] -= a * v;
      if (is_finite(rhi[i])) rhi[i] -= a * v;
      --row_len[i];
    }
    col_removed_[j] = 1;
    fixed_val_[j] = v;
    lb[j] = ub[j] = v;
  };
  auto remove_row = [&](Index i) {
    row_removed[i] = 1;
    for (const auto& [j, a] : row_a[i])
      if (!col_removed_[j]) --col_len[j];
  };
  // f(v) = c v + 1/2 q v^2 of a column on its own.
  auto f1 = [](double cj, double q, double v) { return cj * v + 0.5 * q * v * v; };

  long long changes = 1;
  for (int round = 0; round < kMaxRounds && changes > 0 && !infeasible; ++round) {
    changes = 0;
    ++stats_.rounds;
    // Columns: fixed, empty, dual fixing.
    for (Index j = 0; j < n && !infeasible; ++j) {
      if (col_removed_[j]) continue;
      if (lb[j] > ub[j]) {
        if (lb[j] > ub[j] + tol_of(ub[j])) {
          infeasible = true;
          break;
        }
        lb[j] = ub[j] = p.is_int[j] ? std::round(0.5 * (lb[j] + ub[j])) : 0.5 * (lb[j] + ub[j]);
      }
      if (lb[j] == ub[j]) {
        fix_col(j, lb[j]);
        ++stats_.fixed_cols, ++changes;
        continue;
      }
      if (col_len[j] == 0 && q_len[j] == 0 && qdiag[j] >= 0.0) {  // (a concave term is reported later)
        const double q = qdiag[j];
        double v = 0.0;
        if (q > 0.0) {
          v = std::clamp(-c[j] / q, lb[j], ub[j]);
          if (p.is_int[j]) {
            const double a = std::max(lb[j], std::floor(v)), b = std::min(ub[j], std::ceil(v));
            v = f1(c[j], q, a) <= f1(c[j], q, b) ? a : b;
          }
        } else if (c[j] > 0.0) {
          v = lb[j];
        } else if (c[j] < 0.0) {
          v = ub[j];
        } else {
          v = std::clamp(0.0, lb[j], ub[j]);
        }
        if (!is_finite(v)) continue;  // unbounded along this column if feasible: left to the solver
        fix_col(j, v);
        ++stats_.empty_cols, ++changes;
        continue;
      }
      if (q_len[j] == 0 && qdiag[j] == 0.0 && col_len[j] > 0) {
        // Locks: rows that a decrease (increase) of x_j could violate.
        int down = 0, up = 0;
        for (const auto& [i, a] : col_a[j]) {
          if (row_removed[i]) continue;
          const bool lo = is_finite(rlo[i]), hi = is_finite(rhi[i]);
          if (a > 0.0) down += lo, up += hi;
          else down += hi, up += lo;
        }
        double v = std::numeric_limits<double>::quiet_NaN();
        if (c[j] >= 0.0 && down == 0 && is_finite(lb[j])) v = lb[j];
        else if (c[j] <= 0.0 && up == 0 && is_finite(ub[j])) v = ub[j];
        if (std::isnan(v)) continue;
        fix_col(j, v);
        ++stats_.dual_fixed, ++changes;
      }
    }
    if (infeasible) break;
    // Rows: empty, free, singleton, activity.
    for (Index i = 0; i < m && !infeasible; ++i) {
      if (row_removed[i]) continue;
      if (row_len[i] == 0) {
        if (rlo[i] > tol_of(rlo[i]) || rhi[i] < -tol_of(rhi[i])) {
          infeasible = true;
          break;
        }
        remove_row(i);
        ++stats_.empty_rows, ++changes;
        continue;
      }
      if (!is_finite(rlo[i]) && !is_finite(rhi[i])) {
        remove_row(i);
        ++stats_.redundant_rows, ++changes;
        continue;
      }
      if (row_len[i] == 1) {
        Index j = -1;
        double a = 0.0;
        for (const auto& [k, v] : row_a[i]) {
          if (!col_removed_[k]) {
            j = k, a = v;
            break;
          }
        }
        // Infinite row bounds divide into infinite column bounds of the right sign.
        double l = a > 0.0 ? rlo[i] / a : rhi[i] / a, u = a > 0.0 ? rhi[i] / a : rlo[i] / a;
        if (p.is_int[j]) {
          if (is_finite(l)) l = std::ceil(l - kFeasTol);
          if (is_finite(u)) u = std::floor(u + kFeasTol);
        }
        if (l > lb[j]) lb[j] = l, lb_row_[j] = i;
        if (u < ub[j]) ub[j] = u, ub_row_[j] = i;
        if (lb[j] > ub[j] + tol_of(ub[j])) {
          infeasible = true;
          break;
        }
        singletons_.push_back({i, j, a});
        remove_row(i);
        ++stats_.singleton_rows, ++changes;
        continue;
      }
      // Activity bounds over the current column bounds (counting infinite contributions).
      double amin = 0.0, amax = 0.0;
      int nmin = 0, nmax = 0;
      for (const auto& [j, a] : row_a[i]) {
        if (col_removed_[j]) continue;
        const double lo = a > 0.0 ? lb[j] : ub[j], hi = a > 0.0 ? ub[j] : lb[j];
        if (is_finite(lo)) amin += a * lo;
        else ++nmin;
        if (is_finite(hi)) amax += a * hi;
        else ++nmax;
      }
      if ((nmin == 0 && amin > rhi[i] + tol_of(rhi[i])) || (nmax == 0 && amax < rlo[i] - tol_of(rlo[i]))) {
        infeasible = true;
        break;
      }
      if (is_finite(rlo[i]) && nmin == 0 && amin >= rlo[i] - kRedundantTol * (1.0 + std::abs(rlo[i]))) {
        rlo[i] = -kInf;
        ++changes;
      }
      if (is_finite(rhi[i]) && nmax == 0 && amax <= rhi[i] + kRedundantTol * (1.0 + std::abs(rhi[i]))) {
        rhi[i] = kInf;
        ++changes;
      }
      if (!is_finite(rlo[i]) && !is_finite(rhi[i])) {
        remove_row(i);
        ++stats_.redundant_rows;
        continue;
      }
      // MIQP: integer bounds from the activities (the duals are not needed there).
      if (!has_int) continue;
      for (const auto& [j, a] : row_a[i]) {
        if (col_removed_[j] || !p.is_int[j]) continue;
        const double lo = a > 0.0 ? lb[j] : ub[j], hi = a > 0.0 ? ub[j] : lb[j];
        // Residual activities without column j.
        const bool rmin_ok = nmin == 0 || (nmin == 1 && !is_finite(lo));
        const bool rmax_ok = nmax == 0 || (nmax == 1 && !is_finite(hi));
        const double rmin = nmin == 0 ? amin - a * lo : amin, rmax = nmax == 0 ? amax - a * hi : amax;
        double nl = -kInf, nu = kInf;
        if (is_finite(rhi[i]) && rmin_ok) (a > 0.0 ? nu : nl) = (rhi[i] - rmin) / a;
        if (is_finite(rlo[i]) && rmax_ok) (a > 0.0 ? nl : nu) = (rlo[i] - rmax) / a;
        if (is_finite(nl) && std::abs(nl) < 1e9) {
          nl = std::ceil(nl - kFeasTol);
          if (nl > lb[j]) lb[j] = nl, lb_row_[j] = -1, ++stats_.bound_changes, ++changes;
        }
        if (is_finite(nu) && std::abs(nu) < 1e9) {
          nu = std::floor(nu + kFeasTol);
          if (nu < ub[j]) ub[j] = nu, ub_row_[j] = -1, ++stats_.bound_changes, ++changes;
        }
        if (lb[j] > ub[j]) {
          infeasible = true;
          break;
        }
        // The activities are stale after a change: the row's other columns wait for the next round.
        if ((a > 0.0 ? lb[j] : ub[j]) != lo || (a > 0.0 ? ub[j] : lb[j]) != hi) break;
      }
    }
  }
  if (infeasible) {
    log.log(1, "Presolve: the problem is infeasible\n");
    return QpPresolveStatus::Infeasible;
  }

  // Reduced problem.
  col_map_.clear();
  row_map_.clear();
  std::vector<Index> new_col(static_cast<size_t>(n), -1), new_row(static_cast<size_t>(m), -1);
  for (Index j = 0; j < n; ++j)
    if (!col_removed_[j]) new_col[j] = static_cast<Index>(col_map_.size()), col_map_.push_back(j);
  for (Index i = 0; i < m; ++i)
    if (!row_removed[i]) new_row[i] = static_cast<Index>(row_map_.size()), row_map_.push_back(i);
  const Index rn = static_cast<Index>(col_map_.size()), rm = static_cast<Index>(row_map_.size());
  MipProblem& r = reduced_;
  r = MipProblem();
  r.offset = offset;
  r.display_sense = p.display_sense;
  r.cost.resize(static_cast<size_t>(rn));
  r.col_lb.resize(static_cast<size_t>(rn));
  r.col_ub.resize(static_cast<size_t>(rn));
  r.is_int.resize(static_cast<size_t>(rn));
  std::vector<Index> J, I;
  std::vector<double> V;
  for (Index k = 0; k < rn; ++k) {
    const Index j = col_map_[k];
    r.cost[k] = c[j];
    r.col_lb[k] = lb[j];
    r.col_ub[k] = ub[j];
    r.is_int[k] = p.is_int[j];
    for (const auto& [i, a] : col_a[j]) {
      if (new_row[i] < 0) continue;
      J.push_back(k), I.push_back(new_row[i]), V.push_back(a);
    }
  }
  r.A = SparseMatrix::from_triplets(rn, rm, J, I, V);
  r.row_lb.resize(static_cast<size_t>(rm));
  r.row_ub.resize(static_cast<size_t>(rm));
  for (Index k = 0; k < rm; ++k) r.row_lb[k] = rlo[row_map_[k]], r.row_ub[k] = rhi[row_map_[k]];
  J.clear(), I.clear(), V.clear();
  for (Index j = 0; j < p.Q.num_major; ++j) {
    for (Index k = p.Q.begin(j); k < p.Q.end(j); ++k) {
      const Index i = p.Q.index[k];
      if (new_col[i] < 0 || new_col[j] < 0) continue;
      J.push_back(std::min(new_col[i], new_col[j])), I.push_back(std::max(new_col[i], new_col[j])), V.push_back(p.Q.value[k]);
    }
  }
  r.Q = SparseMatrix::from_triplets(rn, rn, J, I, V);
  stats_.rows_after = rm;
  stats_.cols_after = rn;
  log.log(1,
          "Presolve: %d rows, %d columns -> %d rows, %d columns (%lld fixed, %lld empty and %lld dual-fixed columns; "
          "%lld singleton, %lld redundant and %lld empty rows; %lld integer bounds; %d rounds)\n",
          m, n, rm, rn, stats_.fixed_cols, stats_.empty_cols, stats_.dual_fixed, stats_.singleton_rows,
          stats_.redundant_rows, stats_.empty_rows, stats_.bound_changes, stats_.rounds);
  return QpPresolveStatus::Reduced;
}

std::vector<double> QpPresolver::postsolve(const std::vector<double>& x_reduced) const {
  std::vector<double> x = fixed_val_;
  for (size_t k = 0; k < col_map_.size(); ++k) x[col_map_[k]] = x_reduced[k];
  return x;
}

std::vector<double> QpPresolver::crush(const std::vector<double>& x) const {
  std::vector<double> r(col_map_.size());
  for (size_t k = 0; k < col_map_.size(); ++k) r[k] = x[col_map_[k]];
  return r;
}

void QpPresolver::postsolve_duals(const std::vector<double>& x, const std::vector<double>& y_reduced,
                                  std::vector<double>& y, std::vector<double>& d) const {
  const MipProblem& p = orig_;
  const Index m = p.num_rows();
  y.assign(static_cast<size_t>(m), 0.0);
  for (size_t k = 0; k < row_map_.size() && k < y_reduced.size(); ++k) y[row_map_[k]] = y_reduced[k];
  d = qp_reduced_costs(p, x, y);
  if (singletons_.empty()) return;
  // Singleton rows, last removed first: a row whose bound is its column's final bound on the
  // side the reduced cost presses against takes that reduced cost over (y_i = d_j / a).
  const SparseMatrix At = p.A.transpose();
  for (auto it = singletons_.rbegin(); it != singletons_.rend(); ++it) {
    const Index i = it->row, j = it->col;
    const bool take = (d[j] > 0.0 && lb_row_[j] == i) || (d[j] < 0.0 && ub_row_[j] == i);
    if (!take) continue;
    const double t = d[j] / it->a;
    y[i] += t;
    for (Index k = At.begin(i); k < At.end(i); ++k) d[At.index[k]] -= At.value[k] * t;
    d[j] = 0.0;  // exactly (round-off)
  }
}

}  // namespace milp
