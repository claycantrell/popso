#include "milp/heuristics/feaspump.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <random>
#include <unordered_set>

#include "milp/core/interrupt.hpp"
#include "milp/core/timer.hpp"

namespace milp {

namespace {

uint64_t hash_rounding(const MipProblem& p, const std::vector<double>& xr) {
  uint64_t h = 1469598103934665603ULL;
  for (Index j = 0; j < p.num_cols(); ++j) {
    if (!p.is_int[j]) continue;
    h ^= static_cast<uint64_t>(static_cast<int64_t>(xr[j])) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
  }
  return h;
}

// Moves integer variable j one step away from its current rounded value,
// towards the LP value.
void flip(const MipProblem& p, Index j, double lp_value, std::vector<double>& xr) {
  const double step = lp_value > xr[j] ? 1.0 : -1.0;
  double v = xr[j] + step;
  if (v < p.col_lb[j] || v > p.col_ub[j]) v = xr[j] - step;
  if (v >= p.col_lb[j] && v <= p.col_ub[j]) xr[j] = v;
}

}  // namespace

bool feasibility_pump(const MipProblem& p, const DualSimplex& root, const FeasPumpLimits& limits,
                      std::vector<double>& sol, long long* iterations, std::vector<double>* closest) {
  double best_dist = kInf;
  Timer timer;
  std::mt19937 rng(limits.seed + 17);
  std::uniform_real_distribution<double> unif(0.0, 1.0);
  const Index n = p.num_cols();
  long long iters = 0;
  if (iterations) *iterations = 0;

  std::vector<double> x(root.x(), root.x() + n);
  Basis basis;
  if (root.num_rows() == p.num_rows()) basis = root.get_basis();

  Index num_ints = 0;
  double cnorm = 0.0;
  for (Index j = 0; j < n; ++j) {
    num_ints += p.is_int[j];
    cnorm += p.cost[j] * p.cost[j];
  }
  if (num_ints == 0) return false;
  cnorm = std::sqrt(cnorm);
  const double scale = cnorm > 0 ? std::sqrt(static_cast<double>(num_ints)) / cnorm : 0.0;

  std::vector<double> xr(n), prev(n), cost(n);
  std::unordered_set<uint64_t> seen;
  std::vector<std::pair<double, Index>> order;
  double alpha = 1.0;
  DualSimplex lp;
  // Stall: the projections stopped getting closer to integrality (piperout-08: the best distance
  // came at round 16, then the rounds cycled between the same two roundings until the 10 s time
  // limit). Counted in rounds, so the effort does not depend on the clock or the time limit.
  int stall = 0;
  double stall_best = kInf;
  for (int round = 0; round < limits.max_rounds; ++round) {
    if (timer.elapsed() > limits.time_limit || iters > limits.iteration_budget || interrupted()) break;
    if (limits.stall_rounds > 0 && stall >= limits.stall_rounds) break;
    // Round the LP point.
    for (Index j = 0; j < n; ++j) {
      xr[j] = x[j];
      if (p.is_int[j]) xr[j] = std::min(std::max(std::round(x[j]), p.col_lb[j]), p.col_ub[j]);
    }
    if (round > 0) {
      bool same = true;
      for (Index j = 0; j < n && same; ++j) same = !p.is_int[j] || xr[j] == prev[j];
      if (same) {
        // Short cycle: flip the T integer variables with the largest fractional distance.
        order.clear();
        for (Index j = 0; j < n; ++j) {
          if (p.is_int[j]) order.emplace_back(std::abs(x[j] - xr[j]), j);
        }
        const size_t T = std::min(order.size(), static_cast<size_t>(10 + rng() % 21));
        std::partial_sort(order.begin(), order.begin() + T, order.end(), std::greater<>());
        for (size_t k = 0; k < T; ++k) flip(p, order[k].second, x[order[k].second], xr);
      } else if (seen.count(hash_rounding(p, xr))) {
        // Long cycle: random perturbation.
        for (Index j = 0; j < n; ++j) {
          if (!p.is_int[j]) continue;
          const double r = unif(rng) - 0.3;
          if (std::abs(x[j] - xr[j]) + std::max(r, 0.0) > 0.5) flip(p, j, x[j], xr);
        }
      }
    }
    seen.insert(hash_rounding(p, xr));
    prev = xr;
    if (p.max_violation(xr) <= 1e-6) {
      sol = xr;
      if (iterations) *iterations = iters;
      return true;
    }

    // Project: min distance to the rounding (+ fading objective term).
    alpha *= 0.9;
    for (Index j = 0; j < n; ++j) {
      double d = 0.0;
      if (p.is_int[j]) {
        if (xr[j] <= p.col_lb[j]) d = 1.0;
        else if (xr[j] >= p.col_ub[j]) d = -1.0;
      }
      cost[j] = (1.0 - alpha) * d + alpha * scale * p.cost[j];
    }
    lp.load(p.A, cost, p.col_lb, p.col_ub, p.row_lb, p.row_ub);
    if (!basis.empty()) lp.set_basis(basis);
    LpOptions opt;
    opt.iteration_limit = limits.iteration_budget - iters;
    opt.time_limit = limits.time_limit - timer.elapsed();
    const long long it0 = lp.iterations();
    const LpStatus st = lp.solve(opt);
    iters += lp.iterations() - it0;
    // lp.x() copies (and unscales) all n + m values: take it once per round.
    const double* lx = lp.x();
    if (std::getenv("MILP_FP_DEBUG")) {
      double dist = 0.0;
      for (Index j = 0; j < n; ++j) dist += p.is_int[j] ? std::abs(lx[j] - xr[j]) : 0.0;
      std::fprintf(stderr, "FP round %d status %s dist %.4g viol(xr) %.3g iters %lld t %.2f\n", round, to_string(st),
                   dist, p.max_violation(xr), iters, timer.elapsed());
    }
    if (st != LpStatus::Optimal) break;
    double dist = 0.0;
    for (Index j = 0; j < n; ++j) dist += p.is_int[j] ? std::abs(lx[j] - xr[j]) : 0.0;
    // Improvement: at least 1% closer (tiny decreases while cycling do not count).
    if (dist < 0.99 * stall_best) {
      stall_best = dist;
      stall = 0;
    } else {
      ++stall;
    }
    if (closest) {
      if (dist < best_dist) {
        best_dist = dist;
        closest->assign(lx, lx + n);
      }
    }
    basis = lp.get_basis();
    for (Index j = 0; j < n; ++j) x[j] = lx[j];

    bool integral = true;
    for (Index j = 0; j < n && integral; ++j) integral = !p.is_int[j] || std::abs(x[j] - std::round(x[j])) <= 1e-6;
    if (integral) {
      sol = x;
      for (Index j = 0; j < n; ++j) {
        if (p.is_int[j]) sol[j] = std::round(sol[j]);
      }
      if (p.max_violation(sol) <= 1e-6) {
        if (iterations) *iterations = iters;
        return true;
      }
    }
  }
  if (iterations) *iterations = iters;
  return false;
}

}  // namespace milp
