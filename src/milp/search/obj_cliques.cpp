#include "milp/search/obj_cliques.hpp"

#include <algorithm>
#include <cmath>
#include <utility>

namespace milp {

namespace {

constexpr double kBig = 1e12;  // bounds beyond this count as infinite (no pseudo objective)

}  // namespace

ObjectiveCliques::ObjectiveCliques(const MipProblem& p, long long max_entries)
    : cost_(p.cost), is_int_(p.is_int), n_(p.num_cols()) {
  start_.push_back(0);
  const Index m = p.num_rows();
  const SparseMatrix R = p.A.transpose();
  auto binary = [&](Index j) { return p.is_int[j] && p.col_lb[j] == 0.0 && p.col_ub[j] == 1.0; };
  std::vector<std::pair<double, Index>> items;
  for (Index i = 0; i < m && !truncated_; ++i) {
    if (R.end(i) - R.begin(i) < 2) continue;
    for (int side = 0; side < 2; ++side) {
      const double s = side == 0 ? 1.0 : -1.0;
      const double rhs = side == 0 ? p.row_ub[i] : -p.row_lb[i];
      if (!is_finite(rhs) || std::abs(rhs) > kBig) continue;
      // s a x <= rhs as sum w l <= slack over literals l with weights w > 0, the other columns at their minimum.
      double slack = rhs, mag = std::abs(rhs);
      bool ok = true;
      items.clear();
      for (Index k = R.begin(i); k < R.end(i) && ok; ++k) {
        const Index j = R.index[k];
        const double a = s * R.value[k];
        if (binary(j)) {
          if (a > 0) {
            items.emplace_back(a, 2 * j);
          } else {
            items.emplace_back(-a, 2 * j + 1);
            slack -= a;
          }
          mag += std::abs(a);
        } else {
          const double bnd = a > 0 ? p.col_lb[j] : p.col_ub[j];
          if (!is_finite(bnd) || std::abs(bnd) > kBig) ok = false;
          slack -= a * bnd;
          mag += std::abs(a * bnd);
        }
      }
      if (!ok || items.size() < 2) continue;
      // Pairs must violate the side by more than the feasibility tolerance (as the presolve's clique extension).
      const double tol = 1e-6 * (1.0 + std::abs(rhs)) + 1e-9 * mag;
      std::sort(items.begin(), items.end(), [](const auto& x, const auto& y) {
        return x.first > y.first || (x.first == y.first && x.second < y.second);
      });
      size_t len = 1;
      while (len < items.size() && items[0].first + items[len].first > slack + tol) ++len;
      if (len < 2) continue;
      if (static_cast<long long>(lit_.size() + len) > max_entries) {
        truncated_ = true;
        break;
      }
      for (size_t t = 0; t < len; ++t) {
        lit_.push_back(items[t].second);
        w_.push_back(items[t].first);
        list_of_.push_back(static_cast<Index>(slack_.size()));
      }
      slack_.push_back(slack + tol);
      start_.push_back(static_cast<Index>(lit_.size()));
    }
  }
  // Literal -> entries.
  const size_t nl = 2 * static_cast<size_t>(n_);
  mstart_.assign(nl + 1, 0);
  for (Index l : lit_) ++mstart_[l + 1];
  for (size_t l = 0; l < nl; ++l) mstart_[l + 1] += mstart_[l];
  ment_.resize(lit_.size());
  {
    std::vector<Index> fill(mstart_.begin(), mstart_.end() - 1);
    for (Index e = 0; e < static_cast<Index>(lit_.size()); ++e) ment_[fill[lit_[e]]++] = e;
  }
  // Each literal's increase at p's bounds, all of its lists added up: an upper bound at any tighter bounds (fixed
  // partners add nothing; true ones make the literal infeasible, which this propagation may leave to others).
  inc_.resize(nl);
  for (size_t l = 0; l < nl; ++l) inc_[l] = literal_cost(static_cast<Index>(l), false);
  pre_.resize(lit_.size() + start_.size());
  for (Index r = 0; r < lists(); ++r) {
    Index q = start_[r] + r;
    pre_[q] = 0.0;
    for (Index e = start_[r]; e < start_[r + 1]; ++e, ++q) pre_[q + 1] = pre_[q] + inc_[lit_[e]];
  }
  ub0_.assign(nl, 0.0);
  for (size_t l = 0; l < nl; ++l) {
    if (mstart_[l] == mstart_[l + 1]) continue;
    double u = literal_cost(static_cast<Index>(l), true);
    for (Index k = mstart_[l]; k < mstart_[l + 1]; ++k) {
      const Index e = ment_[k], r = list_of_[e];
      const Index cnt = partners(e);
      u += pre_[start_[r] + r + cnt] - (e < start_[r] + cnt ? inc_[l] : 0.0);
    }
    ub0_[l] = u;
    order_.push_back(static_cast<Index>(l));
  }
  std::sort(order_.begin(), order_.end(), [&](Index a, Index b) { return ub0_[a] > ub0_[b] || (ub0_[a] == ub0_[b] && a < b); });
  list_pass_.assign(start_.size(), -1);
}

// Objective increase over the pseudo objective when literal l is set to 1 (own) or forced to 0, for a free binary.
double ObjectiveCliques::literal_cost(Index l, bool own) const {
  const double c = cost_[l >> 1];
  // Literal 2j (x_j = 1) set: x_j = 1, costs c if c > 0; forced to 0: x_j = 0, costs -c if c < 0. Literal 2j+1
  // (x_j = 0) the other way round.
  const bool up = ((l & 1) == 0) == own;  // the column ends at 1
  return up ? std::max(0.0, c) : std::max(0.0, -c);
}

Index ObjectiveCliques::partners(Index e) const {
  // The entries f with w_f > slack - w_e: a prefix of e's list (e itself included when it is in it).
  const Index r = list_of_[e], b = start_[r], end = start_[r + 1];
  const double thr = slack_[r] - w_[e];
  return static_cast<Index>(std::partition_point(w_.begin() + b, w_.begin() + end, [&](double x) { return x > thr; }) -
                            (w_.begin() + b));
}

Index ObjectiveCliques::propagate(double cutoff, std::vector<double>& lb, std::vector<double>& ub,
                                  std::vector<Index>& changed, long long work_limit, int max_passes,
                                  long long* work_out) {
  long long work = 0;
  Index fixed = 0;
  if (work_out) *work_out = 0;
  if (!is_finite(cutoff) || lit_.empty()) return 0;
  const size_t nl = 2 * static_cast<size_t>(n_);
  true_.resize(nl);
  ptrue_.resize(pre_.size());
  if (stamp_.size() != static_cast<size_t>(n_)) stamp_.assign(n_, 0), stamp_val_.assign(n_, 0);
  auto free_bin = [&](Index j) { return is_int_[j] && lb[j] == 0.0 && ub[j] == 1.0; };
  for (int pass = 0; pass < max_passes && work < work_limit; ++pass) {
    // Pseudo objective at the pass's bounds; every fixing of the pass is derived from these bounds (each is
    // implied by them on its own, so all of them together hold too).
    double P = 0.0, mag = std::abs(cutoff);
    for (Index j = 0; j < n_; ++j) {
      const double c = cost_[j];
      if (c == 0.0) continue;
      const double v = c > 0 ? lb[j] : ub[j];
      if (!(std::abs(v) < kBig)) return fixed;  // no finite pseudo objective
      P += c * v;
      mag += std::abs(c * v);
    }
    work += n_;
    // Rounding errors in P and the increases: fixings need a margin beyond them.
    const double eps = 1e-9 * mag + 1e-7 * std::max(1.0, std::abs(cutoff));
    const double room = cutoff + eps - P;  // a literal is fixed to 0 when its increase exceeds room
    if (room < 0.0) break;  // no improving point within the bounds: left to the LP and the search
    // Candidates: literals whose increase at p's bounds exceeds the room (order_ by that bound, descending).
    size_t ncand = 0;
    while (ncand < order_.size() && ub0_[order_[ncand]] > room) ++ncand;
    if (ncand == 0) break;
    // Objective increase of each literal when forced to 0 (its column moved off its cheaper bound), and whether it
    // is true (the column fixed at the literal's value).
    for (size_t l = 0; l < nl; ++l) {
      const Index j = static_cast<Index>(l >> 1);
      inc_[l] = 0.0;
      true_[l] = 0;
      if (!is_int_[j]) continue;
      if (lb[j] == 0.0 && ub[j] == 1.0) {
        inc_[l] = literal_cost(static_cast<Index>(l), false);
      } else if (lb[j] == ub[j]) {
        true_[l] = lb[j] == ((l & 1) ? 0.0 : 1.0);
      }
    }
    work += n_;
    // Prefix sums of a list (increase, true literals) at pre_ / ptrue_[start_[r] + r + k], built when first needed.
    auto list_prefix = [&](Index r) {
      if (list_pass_[r] == pass_id_) return;
      list_pass_[r] = pass_id_;
      Index q = start_[r] + r;
      pre_[q] = 0.0;
      ptrue_[q] = 0;
      for (Index e = start_[r]; e < start_[r + 1]; ++e, ++q) {
        pre_[q + 1] = pre_[q] + inc_[lit_[e]];
        ptrue_[q + 1] = ptrue_[q] + true_[lit_[e]];
      }
      work += start_[r + 1] - start_[r];
    };
    ++pass_id_;
    Index pass_fixed = 0;
    bool grew = false;
    for (size_t t = 0; t < ncand && work < work_limit; ++t) {
      const Index l = order_[t], v = l >> 1;
      if (!free_bin(v)) continue;
      // Literal l: x_v = 1 (l even) or x_v = 0 (l odd); setting it costs own beyond P.
      const double own = literal_cost(l, true);
      bool fix = own > room;
      double best = 0.0, sum = 0.0;
      for (Index k = mstart_[l]; k < mstart_[l + 1] && !fix; ++k) {
        const Index e = ment_[k], r = list_of_[e], b = start_[r];
        list_prefix(r);
        const Index cnt = partners(e);
        double s = pre_[b + r + cnt];
        int tr = ptrue_[b + r + cnt];
        if (e < b + cnt) s -= inc_[l], tr -= true_[l];
        if (tr > 0) fix = true;  // a partner is true: l = 1 is infeasible
        best = std::max(best, s);
        sum += s;
      }
      // One list alone is a lower bound on the increase, all of them added up an upper bound (a column forced by
      // two lists counts twice); in between, the exact union of the partners.
      if (!fix && own + best <= room && own + sum > room) {
        ++unions_;
        if (++stamp_id_ == 0) std::fill(stamp_.begin(), stamp_.end(), 0), stamp_id_ = 1;
        double tot = own;
        for (Index k = mstart_[l]; k < mstart_[l + 1] && !fix && tot <= room && work < work_limit; ++k) {
          const Index e = ment_[k], r = list_of_[e];
          const double thr = slack_[r] - w_[e];
          for (Index f = start_[r]; f < start_[r + 1] && w_[f] > thr; ++f) {
            ++work;
            if (f == e) continue;
            const Index m = lit_[f], u = m >> 1;
            const char val = static_cast<char>(m & 1);  // which literal of u is forced to 0
            if (u == v) continue;
            if (stamp_[u] == stamp_id_) {
              if (stamp_val_[u] != val) fix = true;  // both values of u excluded: l = 1 is infeasible
              continue;
            }
            stamp_[u] = stamp_id_;
            stamp_val_[u] = val;
            if (true_[m]) fix = true;
            tot += inc_[m];
          }
        }
        fix = fix || tot > room;
      }
      if (!fix && own + best <= room) continue;
      // l = 0 in every improving point: x_v = 0 (l even) or x_v = 1 (l odd).
      const double c = cost_[v];
      if ((l & 1) == 0) {
        ub[v] = 0.0;
        grew = grew || c < 0.0;
      } else {
        lb[v] = 1.0;
        grew = grew || c > 0.0;
      }
      changed.push_back(v);
      ++pass_fixed;
    }
    fixed += pass_fixed;
    if (!grew) break;
  }
  if (work_out) *work_out = work;
  return fixed;
}

}  // namespace milp
