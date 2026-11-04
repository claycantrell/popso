// Tests for the fix-propagate-repair heuristic: every returned point must be feasible (also for
// each strategy of the portfolios alone), runs must be reproducible for a fixed seed and budget,
// and infeasible problems must terminate without a point.

#include <cmath>
#include <random>
#include <string>

#include "milp/core/model.hpp"
#include "milp/heuristics/fpr.hpp"
#include "milp/search/mip_problem.hpp"
#include "test_framework.hpp"

using namespace milp;

namespace {

MipProblem build(const Model& m) {
  MipProblem p;
  bool ok = MipProblem::from_model(m, &p);
  REQUIRE(ok);
  return p;
}

// Random set cover: min c^T x, each row covered by >= 1 of its columns.
MipProblem set_cover(unsigned seed, int rows, int cols) {
  std::mt19937 rng(seed);
  ModelBuilder b;
  for (int j = 0; j < cols; ++j) b.add_col("x" + std::to_string(j), 1.0 + rng() % 20, 0, 1, VarType::Binary);
  for (int i = 0; i < rows; ++i) {
    const Index r = b.add_row("r" + std::to_string(i), 1.0, kInf);
    for (int t = 0; t < 4; ++t) b.add_coef(r, static_cast<Index>(rng() % cols), 1.0);
  }
  return build(b.build());
}

// Random generalized assignment: each job to exactly one machine (equality cliques), machine
// capacities (knapsack rows). Feasible by construction (a random assignment fits).
MipProblem gap(unsigned seed, int jobs, int machines) {
  std::mt19937 rng(seed);
  ModelBuilder b;
  std::vector<int> assign(jobs);
  std::vector<std::vector<int>> w(jobs, std::vector<int>(machines));
  std::vector<double> load(machines, 0.0);
  for (int j = 0; j < jobs; ++j) {
    assign[j] = static_cast<int>(rng() % machines);
    for (int k = 0; k < machines; ++k) w[j][k] = 1 + static_cast<int>(rng() % 9);
    load[assign[j]] += w[j][assign[j]];
  }
  for (int j = 0; j < jobs; ++j)
    for (int k = 0; k < machines; ++k)
      b.add_col("x" + std::to_string(j) + "_" + std::to_string(k), static_cast<double>(rng() % 30), 0, 1,
                VarType::Binary);
  for (int j = 0; j < jobs; ++j) {
    const Index r = b.add_row("job" + std::to_string(j), 1.0, 1.0);
    for (int k = 0; k < machines; ++k) b.add_coef(r, j * machines + k, 1.0);
  }
  for (int k = 0; k < machines; ++k) {
    const Index r = b.add_row("cap" + std::to_string(k), -kInf, load[k] + rng() % 3);
    for (int j = 0; j < jobs; ++j) b.add_coef(r, j * machines + k, w[j][k]);
  }
  return build(b.build());
}

// Random mixed problem: general integers, binaries and bounded continuous variables with
// equality, ranged and one-sided rows, built around a random point (so feasible).
MipProblem mixed(unsigned seed) {
  std::mt19937 rng(seed);
  ModelBuilder b;
  const int n = 12 + static_cast<int>(rng() % 10);
  std::vector<double> x0(n);
  for (int j = 0; j < n; ++j) {
    const int kind = static_cast<int>(rng() % 3);
    const double c = static_cast<double>(static_cast<int>(rng() % 21) - 10);
    if (kind == 0) b.add_col("b" + std::to_string(j), c, 0, 1, VarType::Binary), x0[j] = rng() % 2;
    else if (kind == 1) b.add_col("i" + std::to_string(j), c, -3, 5, VarType::Integer), x0[j] = -3.0 + rng() % 9;
    else b.add_col("c" + std::to_string(j), c, 0, 4), x0[j] = (rng() % 400) / 100.0;
  }
  const int m = 6 + static_cast<int>(rng() % 6);
  for (int i = 0; i < m; ++i) {
    std::vector<std::pair<Index, double>> terms;
    double act = 0.0;
    for (int t = 0; t < 4; ++t) {
      const Index j = static_cast<Index>(rng() % n);
      const double a = static_cast<double>(static_cast<int>(rng() % 11) - 5);
      if (a == 0.0) continue;
      terms.emplace_back(j, a);
      act += a * x0[j];
    }
    const int type = static_cast<int>(rng() % 3);
    const double lo = type == 1 ? -kInf : act - (type == 2 ? 0.0 : rng() % 3);
    const double up = type == 0 ? kInf : act + (type == 2 ? 0.0 : rng() % 3);
    const Index r = b.add_row("r" + std::to_string(i), lo, up);
    for (auto& [j, a] : terms) b.add_coef(r, j, a);
  }
  return build(b.build());
}

std::vector<FprStrategy> all_strategies() {
  std::vector<FprStrategy> s = fpr_lp_free_portfolio();
  s.insert(s.end(), fpr_lp_portfolio().begin(), fpr_lp_portfolio().end());
  return s;
}

}  // namespace

TEST(fpr_set_cover_feasible) {
  for (unsigned seed = 1; seed <= 20; ++seed) {
    MipProblem p = set_cover(seed, 60, 40);
    FprLimits lim;
    lim.effort = 5000000;
    lim.seed = seed;
    std::vector<double> x;
    FprStats st;
    CHECK(fix_propagate_repair(p, lim, x, &st));
    CHECK(p.max_violation(x) <= 1e-6);
    CHECK(st.strategy != nullptr);
  }
}

TEST(fpr_gap_feasible) {
  int found = 0;
  for (unsigned seed = 1; seed <= 20; ++seed) {
    MipProblem p = gap(seed, 30, 5);
    FprLimits lim;
    lim.effort = 20000000;
    lim.seed = seed;
    std::vector<double> x;
    if (fix_propagate_repair(p, lim, x)) {
      ++found;
      CHECK(p.max_violation(x) <= 1e-6);
    }
  }
  CHECK(found >= 15);  // tight capacities make a few hard; most should be found
}

TEST(fpr_mixed_every_point_feasible) {
  int found = 0;
  for (unsigned seed = 1; seed <= 60; ++seed) {
    MipProblem p = mixed(seed);
    FprLimits lim;
    lim.effort = 2000000;
    lim.seed = seed;
    std::vector<double> x;
    if (fix_propagate_repair(p, lim, x)) {
      ++found;
      CHECK(p.max_violation(x) <= 1e-6);
    }
  }
  CHECK(found >= 40);  // all are feasible; FPR is incomplete but should find most
}

TEST(fpr_each_strategy_feasible) {
  // Each strategy alone (incl. the LP-guided ones with an arbitrary reference point) on a mix of
  // problems: whatever it returns must be feasible.
  for (const FprStrategy& s : all_strategies()) {
    for (unsigned seed = 1; seed <= 12; ++seed) {
      MipProblem p = seed % 3 == 0 ? set_cover(seed, 40, 30) : seed % 3 == 1 ? gap(seed, 15, 4) : mixed(seed);
      std::vector<double> ref(p.num_cols());
      for (Index j = 0; j < p.num_cols(); ++j) ref[j] = 0.5 * (std::max(p.col_lb[j], -3.0) + std::min(p.col_ub[j], 3.0));
      FprLimits lim;
      lim.effort = 2000000;
      lim.seed = seed;
      lim.lp_x = &ref;
      lim.strategies = {s};
      std::vector<double> x;
      if (fix_propagate_repair(p, lim, x)) CHECK(p.max_violation(x) <= 1e-6);
    }
  }
}

TEST(fpr_deterministic) {
  for (unsigned seed : {3u, 11u}) {
    MipProblem p = gap(seed, 40, 6);
    auto run = [&](std::vector<double>& x, FprStats& st) {
      FprLimits lim;
      lim.effort = 3000000;
      lim.seed = seed;
      return fix_propagate_repair(p, lim, x, &st);
    };
    std::vector<double> a, b;
    FprStats sa, sb;
    const bool ra = run(a, sa), rb = run(b, sb);
    CHECK(ra == rb);
    CHECK(a == b);
    CHECK(sa.effort == sb.effort);
    CHECK(sa.nodes == sb.nodes);
    CHECK(sa.repairs == sb.repairs);
    CHECK(sa.runs == sb.runs);
  }
}

TEST(fpr_infeasible_terminates) {
  // x + y >= 3 with binaries: root propagation proves infeasibility.
  ModelBuilder b;
  b.add_col("x", 1, 0, 1, VarType::Binary);
  b.add_col("y", 1, 0, 1, VarType::Binary);
  const Index r = b.add_row("r", 3, kInf);
  b.add_coef(r, 0, 1);
  b.add_coef(r, 1, 1);
  MipProblem p = build(b.build());
  FprLimits lim;
  lim.effort = 100000;
  std::vector<double> x;
  FprStats st;
  CHECK(!fix_propagate_repair(p, lim, x, &st));
  CHECK(st.infeasible);

  // Parity: 2 x1 + 2 x2 + 2 x3 = 3 has no integer point, but propagation alone does not see it;
  // the search must stop on its limits.
  ModelBuilder c;
  for (int j = 0; j < 3; ++j) c.add_col("x" + std::to_string(j), 0, 0, 1, VarType::Binary);
  const Index q = c.add_row("q", 3, 3);
  for (int j = 0; j < 3; ++j) c.add_coef(q, j, 2);
  MipProblem p2 = build(c.build());
  CHECK(!fix_propagate_repair(p2, lim, x, &st));
}

TEST(fpr_continuous_completion) {
  // Fixed-charge flow: y_k <= 10 z_k, sum y_k >= 25, costs on z and y; the continuous part comes
  // from the leaf LP.
  ModelBuilder b;
  for (int k = 0; k < 4; ++k) b.add_col("z" + std::to_string(k), 5.0 + k, 0, 1, VarType::Binary);
  for (int k = 0; k < 4; ++k) b.add_col("y" + std::to_string(k), 1.0, 0, kInf);
  for (int k = 0; k < 4; ++k) {
    const Index r = b.add_row("link" + std::to_string(k), -kInf, 0);
    b.add_coef(r, 4 + k, 1);
    b.add_coef(r, k, -10);
  }
  const Index d = b.add_row("demand", 25, kInf);
  for (int k = 0; k < 4; ++k) b.add_coef(d, 4 + k, 1);
  MipProblem p = build(b.build());
  FprLimits lim;
  lim.effort = 1000000;
  std::vector<double> x;
  FprStats st;
  CHECK(fix_propagate_repair(p, lim, x, &st));
  CHECK(p.max_violation(x) <= 1e-6);
}
