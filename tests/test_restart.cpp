// Restarts: pseudocosts carried through presolve (mapping, capping), repeated root restarts (forced
// by low thresholds) that keep the optimum of a solve without restarts in both modes, with
// reproducible deterministic reruns, and the restart budget.

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <initializer_list>
#include <random>
#include <string>
#include <utility>
#include <vector>

#include "milp/core/log.hpp"
#include "milp/core/model.hpp"
#include "milp/presolve/presolve.hpp"
#include "milp/search/branch_and_bound.hpp"
#include "milp/search/mip_problem.hpp"
#include "milp/solver/solver.hpp"
#include "test_framework.hpp"

using namespace milp;

namespace {

void set_env_value(const char* name, const char* value) {
#ifdef _WIN32
  _putenv_s(name, value ? value : "");
#else
  if (value) setenv(name, value, 1);
  else unsetenv(name);
#endif
}

// Sets restart knobs for the lifetime of the object (and restores the defaults after).
struct RestartEnv {
  std::vector<const char*> names;
  RestartEnv(std::initializer_list<std::pair<const char*, const char*>> kv) {
    for (const auto& [k, v] : kv) {
      set_env_value(k, v);
      names.push_back(k);
    }
  }
  ~RestartEnv() {
    for (const char* k : names) set_env_value(k, nullptr);
  }
};

// Snapshot entry values that identify the column: sum 10 j + d + 0.5, count j % 7 + 1 (direction d).
PseudocostSnapshot labelled_snapshot(Index n) {
  PseudocostSnapshot s;
  for (Index j = 0; j < n; ++j) {
    for (int d = 0; d < 2; ++d) {
      s.sum.push_back(10.0 * j + d + 0.5);
      s.cnt.push_back(static_cast<int>(j % 7) + 1);
      s.inf_sum.push_back(100.0 * j + d);
      s.inf_cnt.push_back(static_cast<int>(j % 5) + 1);
    }
  }
  return s;
}

}  // namespace

// Pseudocosts through presolve: kept columns keep their entries, removed (fixed) columns disappear,
// a merged parallel column's surviving column keeps its own entry.
TEST(restart_pseudocosts_mapped_through_presolve) {
  ModelBuilder b;
  // x0, x1: parallel (same column and cost): merged into one of them. x2: fixed by its bounds.
  // x3, x4: kept. c5: continuous, kept.
  const Index x0 = b.add_col("x0", 2.0, 0.0, 1.0, VarType::Integer);
  const Index x1 = b.add_col("x1", 2.0, 0.0, 1.0, VarType::Integer);
  const Index x2 = b.add_col("x2", 1.0, 1.0, 1.0, VarType::Integer);
  const Index x3 = b.add_col("x3", 3.0, 0.0, 4.0, VarType::Integer);
  const Index x4 = b.add_col("x4", -1.0, 0.0, 3.0, VarType::Integer);
  const Index c5 = b.add_col("c5", 1.5, 0.0, 10.0);
  const Index r0 = b.add_row("r0", 2.0, kInf);
  for (auto [j, v] : std::vector<std::pair<Index, double>>{{x0, 1}, {x1, 1}, {x2, 1}, {x3, 2}, {x4, 1}, {c5, 1}})
    b.add_coef(r0, j, v);
  const Index r1 = b.add_row("r1", -kInf, 7.0);
  for (auto [j, v] : std::vector<std::pair<Index, double>>{{x0, 2}, {x1, 2}, {x3, 1}, {x4, 3}, {c5, -1}})
    b.add_coef(r1, j, v);
  const Index r2 = b.add_row("r2", -kInf, 5.0);
  for (auto [j, v] : std::vector<std::pair<Index, double>>{{x3, 1}, {x4, 1}, {c5, 1}}) b.add_coef(r2, j, v);
  MipProblem p;
  REQUIRE(MipProblem::from_model(b.build(), &p));
  PresolveOptions o;
  o.probing = false;
  o.dominated_cols = false;
  o.dual_fixing = false;
  Presolver pre;
  REQUIRE(pre.run(p, o, Logger(0)) == PresolveStatus::Reduced);
  const std::vector<Index>& map = pre.col_map();
  REQUIRE(map.size() < static_cast<size_t>(p.num_cols()));  // x2 fixed, one of x0 / x1 merged
  bool has_x2 = false, has_x0 = false, has_x1 = false;
  for (Index j : map) has_x2 |= j == x2, has_x0 |= j == x0, has_x1 |= j == x1;
  CHECK(!has_x2);
  CHECK(has_x0 != has_x1);
  const PseudocostSnapshot s = labelled_snapshot(p.num_cols());
  const PseudocostSnapshot t = s.mapped(map);
  REQUIRE(t.sum.size() == 2 * map.size());
  REQUIRE(t.cnt.size() == 2 * map.size() && t.inf_sum.size() == 2 * map.size() && t.inf_cnt.size() == 2 * map.size());
  for (size_t k = 0; k < map.size(); ++k) {
    const size_t j = static_cast<size_t>(map[k]);
    for (size_t d = 0; d < 2; ++d) {
      CHECK_EQ(t.sum[2 * k + d], s.sum[2 * j + d]);
      CHECK_EQ(t.cnt[2 * k + d], s.cnt[2 * j + d]);
      CHECK_EQ(t.inf_sum[2 * k + d], s.inf_sum[2 * j + d]);
      CHECK_EQ(t.inf_cnt[2 * k + d], s.inf_cnt[2 * j + d]);
    }
  }
  CHECK(PseudocostSnapshot{}.mapped(map).empty());
  // A start crushed the same way lands on the same columns (the mapping restarts use for both).
  std::vector<double> x(static_cast<size_t>(p.num_cols()));
  for (Index j = 0; j < p.num_cols(); ++j) x[j] = static_cast<double>(j);
  const std::vector<double> xc = pre.crush(x);
  for (size_t k = 0; k < map.size(); ++k) {
    if (map[k] == x3 || map[k] == x4 || map[k] == c5) CHECK_EQ(xc[k], static_cast<double>(map[k]));
  }
}

// Capping keeps the averages and limits the counts to the reliability threshold.
TEST(restart_pseudocosts_capped) {
  const PseudocostSnapshot s = labelled_snapshot(9);
  const PseudocostSnapshot c = s.capped(4);
  REQUIRE(c.cnt.size() == s.cnt.size());
  for (size_t k = 0; k < s.cnt.size(); ++k) {
    CHECK(c.cnt[k] == std::min(s.cnt[k], 4));
    CHECK_NEAR(c.sum[k] / c.cnt[k], s.sum[k] / s.cnt[k], 1e-12);
    CHECK(c.inf_cnt[k] == std::min(s.inf_cnt[k], 4));
    CHECK_NEAR(c.inf_sum[k] / c.inf_cnt[k], s.inf_sum[k] / s.inf_cnt[k], 1e-12);
  }
}

namespace {

// Multi-dimensional knapsack (maximize) with values correlated to the weights: root reduced-cost
// fixing with a good incumbent fixes many columns, at every root.
Model restart_mkp(unsigned seed, int n, int m) {
  std::mt19937 rng(seed);
  ModelBuilder b;
  b.model().sense = ObjSense::Maximize;
  std::vector<std::vector<int>> a(m, std::vector<int>(n));
  std::vector<double> sum(m, 0.0);
  for (int i = 0; i < m; ++i) {
    for (int j = 0; j < n; ++j) sum[i] += a[i][j] = 10 + static_cast<int>(rng() % 91);
  }
  for (int j = 0; j < n; ++j) {
    double v = static_cast<int>(rng() % 41);
    for (int i = 0; i < m; ++i) v += a[i][j] / static_cast<double>(m);
    b.add_col("x" + std::to_string(j), std::round(v), 0.0, 1.0, VarType::Binary);
  }
  for (int i = 0; i < m; ++i) {
    const Index row = b.add_row("c" + std::to_string(i), -kInf, std::floor(0.4 * sum[i]));
    for (int j = 0; j < n; ++j) b.add_coef(row, j, a[i][j]);
  }
  return b.build();
}

SolveResult solve_restart_model(const Model& m, int threads, bool det, unsigned seed) {
  SolverParams sp;
  sp.verbosity = std::getenv("MILP_TEST_VERBOSE") ? 1 : 0;
  sp.threads = threads;
  sp.deterministic = det;
  sp.seed = seed;
  sp.mip_rel_gap = 0.0;
  sp.mip_abs_gap = 1e-9;
  sp.node_limit = 500000;
  return solve(m, sp);
}

}  // namespace

// Root restarts forced at every root (any newly fixed column, unlimited budget): several restarts
// in a row, and the optimum of a solve without restarts, in both modes; deterministic reruns are
// identical. Pseudocosts carried as history (seed 3), as estimates (seed 5) and not (seed 8).
TEST(restart_repeated_root_restarts_keep_optimum) {
  int multi = 0;
  for (unsigned seed : {3u, 5u, 8u}) {
    const char* pc_mode = seed == 3u ? "1" : seed == 5u ? "2" : "0";
    const Model m = restart_mkp(seed, 120, 3);
    SolveResult ref;
    {
      RestartEnv env({{"MILP_MAX_RESTARTS", "0"}, {"MILP_MAX_TREE_RESTARTS", "0"}});
      ref = solve_restart_model(m, 1, false, 1);
    }
    REQUIRE(ref.status == Status::Optimal);
    CHECK_EQ(ref.stats.restarts, 0);
    RestartEnv env({{"MILP_RESTART_FRAC", "0.001"},
                    {"MILP_RESTART_FRAC2", "0.001"},
                    {"MILP_RESTART_MIN", "1"},
                    {"MILP_RESTART_MINRED", "0"},
                    {"MILP_RESTART_WORK", "1e9"},
                    {"MILP_RESTART_PC", pc_mode}});
    const SolveResult a = solve_restart_model(m, 1, false, 1);
    CHECK(a.status == Status::Optimal);
    CHECK_NEAR(a.objective, ref.objective, 1e-6);
    CHECK(a.stats.restarts <= 4);
    const SolveResult d1 = solve_restart_model(m, 2, true, 1);
    const SolveResult d2 = solve_restart_model(m, 2, true, 1);
    CHECK(d1.status == Status::Optimal);
    CHECK_NEAR(d1.objective, ref.objective, 1e-6);
    CHECK(d1.status == d2.status);
    CHECK_EQ(d1.stats.restarts, d2.stats.restarts);
    CHECK_EQ(d1.stats.nodes, d2.stats.nodes);
    CHECK_EQ(d1.stats.lp_iterations, d2.stats.lp_iterations);
    CHECK_EQ(d1.objective, d2.objective);
    multi += a.stats.restarts >= 2 && d1.stats.restarts >= 2;
    if (std::getenv("MILP_TEST_VERBOSE")) {
      std::printf("  mkp seed %u: restarts %d (deterministic %d), nodes %lld / %lld (without restarts %lld)\n", seed,
                  a.stats.restarts, d1.stats.restarts, a.stats.nodes, d1.stats.nodes, ref.stats.nodes);
    }
  }
  CHECK(multi >= 1);  // repeated restarts happened
}

// The budget of the restarted roots: none left (MILP_RESTART_WORK=0) allows the first restart only.
TEST(restart_budget_limits_repeated_restarts) {
  for (unsigned seed : {3u, 5u, 8u}) {
    const Model m = restart_mkp(seed, 120, 3);
    RestartEnv env({{"MILP_RESTART_FRAC", "0.001"},
                    {"MILP_RESTART_FRAC2", "0.001"},
                    {"MILP_RESTART_MIN", "1"},
                    {"MILP_RESTART_MINRED", "0"},
                    {"MILP_RESTART_WORK", "0"}});
    for (bool det : {false, true}) {
      const SolveResult r = solve_restart_model(m, det ? 2 : 1, det, 1);
      CHECK(r.status == Status::Optimal);
      CHECK(r.stats.restarts <= 1);
    }
  }
}

// The restart rule shared by the end of the root and root strong branching (root_restart_decision):
// the share and minimum count of newly fixed integers, and from the second restart on the refusals
// (the last presolve removed too little; the budget, in work units deterministically, else seconds).
TEST(restart_decision_rule) {
  SearchOptions o;
  o.restart_frac = 0.1;
  o.restart_min_cols = 10;
  o.restart_min_removed = 0.05;
  using D = RestartDecision;
  // First restart: no budget, no removal check.
  CHECK(root_restart_decision(o, false, 200, 19, 1e12, 1e12).outcome == D::TooFew);
  CHECK_EQ(root_restart_decision(o, false, 200, 19, 0, 0).need, 20);
  CHECK(root_restart_decision(o, false, 200, 20, 1e12, 1e12).outcome == D::Restart);
  CHECK(root_restart_decision(o, true, 50, 9, 0, 0).outcome == D::TooFew);  // at least restart_min_cols
  CHECK(root_restart_decision(o, true, 50, 10, 0, 0).outcome == D::Restart);
  CHECK(root_restart_decision(o, true, 0, 0, 0, 0).outcome == D::TooFew);
  // Chained restarts.
  o.restarts_done = 1;
  o.restart_work_left = 100.0;
  o.restart_time_left = 5.0;
  o.last_restart_removed = 0.04;
  D d = root_restart_decision(o, true, 200, 50, 10, 0.1);
  CHECK(d.outcome == D::Refused);
  REQUIRE(d.refused != nullptr);
  CHECK(std::string(d.refused).find("removed") != std::string::npos);
  o.last_restart_removed = 0.05;
  CHECK(root_restart_decision(o, true, 200, 50, 50, 1e9).outcome == D::Restart);  // 2 x 50 <= 100
  d = root_restart_decision(o, true, 200, 50, 51, 0.0);
  CHECK(d.outcome == D::Refused);
  REQUIRE(d.refused != nullptr);
  CHECK(std::string(d.refused).find("budget") != std::string::npos);
  CHECK(root_restart_decision(o, false, 200, 50, 1e12, 2.5).outcome == D::Restart);  // seconds: 2 x 2.5 <= 5
  CHECK(root_restart_decision(o, false, 200, 50, 0.0, 2.6).outcome == D::Refused);
  CHECK(root_restart_decision(o, false, 200, 19, 0.0, 0.0).outcome == D::TooFew);  // too few comes first
}

// Root strong branching skipped in roots that restart anyway (MILP_RESTART_SKIP_SB=1), with restarts
// forced at every root: the optimum of a solve without restarts, and identical deterministic reruns.
TEST(restart_skip_root_strong_branching_keeps_optimum) {
  for (unsigned seed : {3u, 5u, 8u}) {
    const Model m = restart_mkp(seed, 120, 3);
    SolveResult ref;
    {
      RestartEnv env({{"MILP_MAX_RESTARTS", "0"}, {"MILP_MAX_TREE_RESTARTS", "0"}});
      ref = solve_restart_model(m, 1, false, 1);
    }
    REQUIRE(ref.status == Status::Optimal);
    RestartEnv env({{"MILP_RESTART_FRAC", "0.001"},
                    {"MILP_RESTART_FRAC2", "0.001"},
                    {"MILP_RESTART_MIN", "1"},
                    {"MILP_RESTART_MINRED", "0"},
                    {"MILP_RESTART_WORK", "1e9"},
                    {"MILP_RESTART_SKIP_SB", "1"},
                    {"MILP_RESTART_SKIP_SB_WORK", "0"}});
    const SolveResult a = solve_restart_model(m, 1, false, 1);
    CHECK(a.status == Status::Optimal);
    CHECK_NEAR(a.objective, ref.objective, 1e-6);
    const SolveResult d1 = solve_restart_model(m, 2, true, 1);
    const SolveResult d2 = solve_restart_model(m, 2, true, 1);
    CHECK(d1.status == Status::Optimal);
    CHECK_NEAR(d1.objective, ref.objective, 1e-6);
    CHECK_EQ(d1.stats.restarts, d2.stats.restarts);
    CHECK_EQ(d1.stats.nodes, d2.stats.nodes);
    CHECK_EQ(d1.stats.lp_iterations, d2.stats.lp_iterations);
    CHECK_EQ(d1.objective, d2.objective);
    if (std::getenv("MILP_TEST_VERBOSE")) {
      std::printf("  mkp seed %u, root strong branching skipped: restarts %d (deterministic %d), nodes %lld / %lld\n",
                  seed, a.stats.restarts, d1.stats.restarts, a.stats.nodes, d1.stats.nodes);
    }
  }
}
