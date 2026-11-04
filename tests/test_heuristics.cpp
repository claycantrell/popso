// Tests for the primal heuristics: Feasibility Jump, 1-opt, local branching,
// mutation and propagating dives. Every reported point must be feasible; FJ
// must be reproducible.

#include <cmath>
#include <random>
#include <string>

#include "milp/core/model.hpp"
#include "milp/heuristics/diving.hpp"
#include "milp/heuristics/fj.hpp"
#include "milp/heuristics/lmip.hpp"
#include "milp/heuristics/localsearch.hpp"
#include "milp/heuristics/submip.hpp"
#include "milp/search/mip_problem.hpp"
#include "test_framework.hpp"

using namespace milp;

namespace {

MipProblem make_problem(const Model& m) {
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
  return make_problem(b.build());
}

// Random mixed problem: general integers, binaries and bounded continuous
// variables with equality, ranged and one-sided rows. Feasible by construction
// (rows are built around a random point).
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
  return make_problem(b.build());
}

}  // namespace

TEST(fj_set_cover_feasible_and_improving) {
  for (unsigned seed = 1; seed <= 20; ++seed) {
    MipProblem p = set_cover(seed, 60, 40);
    FjLimits lim;
    lim.effort = 2000000;
    lim.seed = seed;
    double last = kInf;
    int sols = 0;
    bool ok = feasibility_jump(p, lim, [&](const std::vector<double>& x, double obj) {
      CHECK(p.max_violation(x) <= 1e-6);
      CHECK(obj < last);
      CHECK_NEAR(obj, p.objective(x) - p.offset, 1e-9);
      last = obj;
      ++sols;
    });
    CHECK(ok);
    CHECK(sols >= 1);
  }
}

TEST(fj_mixed_feasible) {
  int found = 0;
  for (unsigned seed = 1; seed <= 60; ++seed) {
    MipProblem p = mixed(seed);
    FjLimits lim;
    lim.effort = 1000000;
    lim.seed = seed;
    found += feasibility_jump(p, lim, [&](const std::vector<double>& x, double) {
      CHECK(p.max_violation(x) <= 1e-6);
    });
  }
  CHECK(found >= 40);  // all are feasible; FJ is incomplete but should find most
}

TEST(fj_deterministic) {
  MipProblem p = set_cover(7, 200, 120);
  auto run = [&](std::vector<double>& objs) {
    FjLimits lim;
    lim.effort = 3000000;
    lim.seed = 3;
    FjStats st;
    feasibility_jump(p, lim, [&](const std::vector<double>&, double obj) { objs.push_back(obj); }, &st);
    objs.push_back(static_cast<double>(st.moves));
  };
  std::vector<double> a, b;
  run(a);
  run(b);
  CHECK(a == b);
}

TEST(fj_infeasible_terminates) {
  // x + y >= 3 with binaries: infeasible; FJ must stop on its budget.
  ModelBuilder b;
  b.add_col("x", 1, 0, 1, VarType::Binary);
  b.add_col("y", 1, 0, 1, VarType::Binary);
  const Index r = b.add_row("r", 3, kInf);
  b.add_coef(r, 0, 1);
  b.add_coef(r, 1, 1);
  MipProblem p = make_problem(b.build());
  FjLimits lim;
  lim.effort = 100000;
  CHECK(!feasibility_jump(p, lim, [](const std::vector<double>&, double) {}));
}

TEST(one_opt_keeps_feasibility) {
  for (unsigned seed = 1; seed <= 20; ++seed) {
    MipProblem p = set_cover(seed, 30, 30);
    std::vector<double> x(p.num_cols(), 1.0);  // all-ones covers everything
    REQUIRE(p.max_violation(x) <= 1e-9);
    const double before = p.objective(x);
    CHECK(one_opt(p, x));
    CHECK(p.max_violation(x) <= 1e-9);
    CHECK(p.objective(x) < before);
  }
  for (unsigned seed = 1; seed <= 60; ++seed) {
    MipProblem p = mixed(seed);
    FjLimits lim;
    lim.effort = 1000000;
    std::vector<double> x;
    feasibility_jump(p, lim, [&](const std::vector<double>& s, double) { x = s; });
    if (x.empty()) continue;
    const double before = p.objective(x);
    one_opt(p, x);
    CHECK(p.max_violation(x) <= 1e-6);
    CHECK(p.objective(x) <= before + 1e-9);
  }
}

TEST(local_branching_neighbourhood) {
  MipProblem p = set_cover(11, 40, 30);
  std::vector<double> inc(p.num_cols(), 1.0);
  SubMipLimits lim;
  lim.node_limit = 1000;
  std::vector<double> sol;
  REQUIRE(local_branching(p, inc, 5, p.objective(inc), lim, sol));
  CHECK(p.max_violation(sol) <= 1e-6);
  CHECK(p.objective(sol) < p.objective(inc));
  double dist = 0.0;
  for (Index j = 0; j < p.num_cols(); ++j) dist += std::abs(sol[j] - inc[j]);
  CHECK(dist <= 5.0 + 1e-6);
}

TEST(fj_row_of_fixed_columns_terminates) {
  // x fixed at 0 makes x + y >= 2 unsatisfiable (y <= 1): FJ must stop on its effort budget.
  ModelBuilder b;
  b.add_col("x", 1, 0, 0, VarType::Integer);
  b.add_col("y", 1, 0, 1, VarType::Binary);
  const Index r = b.add_row("r", 2, kInf);
  b.add_coef(r, 0, 1);
  b.add_coef(r, 1, 1);
  const Index r2 = b.add_row("only_fixed", 1, kInf);
  b.add_coef(r2, 0, 1);
  MipProblem p = make_problem(b.build());
  FjLimits lim;
  lim.effort = 100000;
  CHECK(!feasibility_jump(p, lim, [](const std::vector<double>&, double) {}));
}

TEST(mutation_neighbourhood) {
  MipProblem p = set_cover(13, 40, 30);
  std::vector<double> inc(p.num_cols(), 1.0);
  for (double frac : {0.3, 0.6}) {
    SubMipStats st;
    SubMipLimits lim;
    lim.node_limit = 1000;
    lim.stats = &st;
    std::vector<double> sol;
    REQUIRE(mutation(p, inc, frac, 5, p.objective(inc), lim, sol));
    CHECK(p.max_violation(sol) <= 1e-6);
    CHECK(p.objective(sol) < p.objective(inc));
    CHECK(st.solves == 1);
    CHECK(st.fixed == static_cast<Index>(frac * p.num_cols()));
  }
  // Nothing (or everything) to fix: rejected without a solve.
  SubMipStats st;
  SubMipLimits lim;
  lim.stats = &st;
  std::vector<double> sol;
  CHECK(!mutation(p, inc, 1.0, 5, kInf, lim, sol));
  CHECK(st.solves == 0);
}

TEST(propagating_dive_feasible_and_restores_bounds) {
  int found = 0;
  for (unsigned seed = 1; seed <= 40; ++seed) {
    MipProblem p = mixed(seed);
    DualSimplex lp;
    lp.load(p.A, p.cost, p.col_lb, p.col_ub, p.row_lb, p.row_ub);
    if (lp.solve(LpOptions{}) != LpStatus::Optimal) continue;
    SimpleRounding rounding(p);
    DomainPropagator prop(p);
    for (DiveRule rule : {DiveRule::Fractional, DiveRule::Coefficient, DiveRule::VectorLength}) {
      REQUIRE(lp.solve(LpOptions{}) == LpStatus::Optimal);  // a dive leaves the LP where it ended
      std::vector<double> sol;
      long long its = 0;
      if (dive(p, lp, rule, rounding, nullptr, DiveLimits{}, sol, &its, &prop)) {
        CHECK(p.max_violation(sol) <= 1e-6);
        ++found;
      }
      for (Index j = 0; j < p.num_cols(); ++j) {
        CHECK(lp.col_lower(j) == p.col_lb[j]);
        CHECK(lp.col_upper(j) == p.col_ub[j]);
      }
    }
  }
  CHECK(found >= 30);
}

TEST(repair_infeasible_point) {
  // All ones except the columns of row 0, which is then uncovered: repair frees
  // that row's neighbourhood and covers it again.
  for (unsigned seed = 1; seed <= 10; ++seed) {
    MipProblem p = set_cover(seed, 40, 30);
    const SparseMatrix rows = p.A.transpose();
    std::vector<double> x(p.num_cols(), 1.0);
    for (Index k = rows.begin(0); k < rows.end(0); ++k) x[rows.index[k]] = 0.0;
    REQUIRE(p.max_violation(x) > 1e-6);
    SubMipStats st;
    SubMipLimits lim;
    lim.node_limit = 1000;
    lim.stats = &st;
    std::vector<double> sol;
    REQUIRE(repair(p, x, kInf, lim, sol));
    CHECK(p.max_violation(sol) <= 1e-6);
    CHECK(st.solves == 1);
    CHECK(st.fixed > 0);
  }
}

TEST(rins_target_fixing) {
  MipProblem p = set_cover(19, 40, 30);
  std::vector<double> inc(p.num_cols(), 1.0), lp(p.num_cols(), 0.5);
  for (Index j = 0; j < 10; ++j) lp[j] = 1.0;  // agrees on 10 of 30 columns
  for (Index j = 10; j < p.num_cols(); ++j) lp[j] = 0.5 + 0.01 * static_cast<double>(j);  // closer later
  SubMipStats st;
  SubMipLimits lim;
  lim.node_limit = 1000;
  lim.stats = &st;
  lim.min_fixed_fraction = 0.0;
  std::vector<double> sol;
  rins(p, lp.data(), inc, p.objective(inc), lim, sol, 0.6);
  CHECK(st.solves == 1);
  CHECK(st.fixed == 18);  // 0.6 * 30: the 10 agreeing plus the 8 closest (columns 22..29)
}

TEST(crossover_fixes_agreeing_columns) {
  // Two feasible covers (all ones, and all ones but the first 8 columns where removable):
  // crossover fixes the columns where they agree and must return a feasible, better point.
  MipProblem p = set_cover(23, 40, 30);
  std::vector<double> a(p.num_cols(), 1.0), b = a;
  for (Index j = 0; j < 8; ++j) {
    b[j] = 0.0;
    if (p.max_violation(b) > 1e-6) b[j] = 1.0;
  }
  Index agree = 0;
  for (Index j = 0; j < p.num_cols(); ++j) agree += a[j] == b[j];
  SubMipStats st;
  SubMipLimits lim;
  lim.node_limit = 1000;
  lim.min_fixed_fraction = 0.0;
  lim.stats = &st;
  std::vector<double> sol;
  const double cut = std::min(p.objective(a), p.objective(b));
  const bool ok = crossover(p, {&a, &b}, cut, lim, sol);
  CHECK(st.solves == 1);
  CHECK(st.fixed == agree);
  if (ok) {
    CHECK(p.max_violation(sol) <= 1e-6);
    CHECK(p.objective(sol) < cut);
    for (Index j = 0; j < p.num_cols(); ++j) {
      if (a[j] == b[j]) CHECK(sol[j] == a[j]);
    }
  }
  // A single point, or identical points (nothing left free), is rejected without a solve.
  SubMipStats st2;
  lim.stats = &st2;
  CHECK(!crossover(p, {&a}, kInf, lim, sol));
  CHECK(!crossover(p, {&a, &a}, kInf, lim, sol));
  CHECK(st2.solves == 0);
}

// Local-MIP style search (heuristics/lmip.hpp): every reported point is feasible, objectives
// strictly decrease, the run is reproducible for a seed and budget.
TEST(lmip_set_cover_feasible_and_improving) {
  for (unsigned seed = 1; seed <= 20; ++seed) {
    MipProblem p = set_cover(seed, 60, 40);
    LocalMipLimits lim;
    lim.effort = 2000000;
    lim.seed = seed;
    double last = kInf;
    int sols = 0;
    bool ok = local_mip(p, lim, [&](const std::vector<double>& x, double obj) {
      CHECK(p.max_violation(x) <= 1e-6);
      CHECK(obj < last);
      CHECK_NEAR(obj, p.objective(x) - p.offset, 1e-9);
      last = obj;
      ++sols;
    });
    CHECK(ok);
    CHECK(sols >= 1);
  }
}

TEST(lmip_mixed_feasible) {
  int found = 0;
  for (unsigned seed = 1; seed <= 60; ++seed) {
    MipProblem p = mixed(seed);
    LocalMipLimits lim;
    lim.effort = 1000000;
    lim.seed = seed;
    found += local_mip(p, lim, [&](const std::vector<double>& x, double obj) {
      CHECK(p.max_violation(x) <= 1e-6);
      CHECK_NEAR(obj, p.objective(x) - p.offset, 1e-6);
    });
  }
  CHECK(found >= 40);
}

TEST(lmip_start_and_cutoff) {
  // From a feasible start with its value as the cutoff, only strictly better points are reported.
  for (unsigned seed = 1; seed <= 10; ++seed) {
    MipProblem p = set_cover(seed, 40, 30);
    std::vector<double> start(p.num_cols(), 1.0);
    const double c0 = p.objective(start) - p.offset;
    LocalMipLimits lim;
    lim.effort = 1000000;
    lim.seed = seed;
    lim.start = &start;
    lim.cutoff = c0;
    bool any = local_mip(p, lim, [&](const std::vector<double>& x, double obj) {
      CHECK(p.max_violation(x) <= 1e-6);
      CHECK(obj < c0);
    });
    CHECK(any);  // all-ones is far from optimal
  }
}

TEST(lmip_deterministic) {
  for (unsigned s : {7u, 8u}) {
    MipProblem p = s == 7 ? set_cover(s, 200, 120) : mixed(s);
    auto run = [&](std::vector<double>& out) {
      LocalMipLimits lim;
      lim.effort = 3000000;
      lim.seed = 3;
      LocalMipStats st;
      local_mip(p, lim, [&](const std::vector<double>& x, double obj) {
        out.push_back(obj);
        out.insert(out.end(), x.begin(), x.end());
      }, &st);
      out.push_back(static_cast<double>(st.steps));
      out.push_back(static_cast<double>(st.effort));
    };
    std::vector<double> a, b;
    run(a);
    run(b);
    CHECK(a == b);
  }
}

TEST(lmip_infeasible_terminates) {
  // x fixed at 0 makes x + y >= 2 unsatisfiable: the search must stop on its effort budget.
  ModelBuilder b;
  b.add_col("x", 1, 0, 0, VarType::Integer);
  b.add_col("y", 1, 0, 1, VarType::Binary);
  const Index r = b.add_row("r", 2, kInf);
  b.add_coef(r, 0, 1);
  b.add_coef(r, 1, 1);
  MipProblem p = make_problem(b.build());
  LocalMipLimits lim;
  lim.effort = 100000;
  CHECK(!local_mip(p, lim, [](const std::vector<double>&, double) {}));
  // With random restarts on stalls it still ends on the effort budget.
  lim.stall_effort = 1000;
  lim.stall_restart = true;
  LocalMipStats st;
  CHECK(!local_mip(p, lim, [](const std::vector<double>&, double) {}, &st));
  CHECK(st.restarts >= 2);
  CHECK(st.effort >= lim.effort && st.effort < lim.effort + 100000);
}

TEST(lmip_stall_restart_mixed_feasible) {
  // Short stalls with restarts: every point reported is still feasible.
  int found = 0;
  for (unsigned seed = 1; seed <= 30; ++seed) {
    MipProblem p = mixed(seed);
    LocalMipLimits lim;
    lim.effort = 500000;
    lim.stall_effort = 5000;
    lim.stall_restart = true;
    lim.seed = seed;
    found += local_mip(p, lim, [&](const std::vector<double>& x, double) { CHECK(p.max_violation(x) <= 1e-6); });
  }
  CHECK(found >= 20);
}

TEST(lp_rens_feasible_and_deterministic) {
  // The root-LP sub-MIP must only report feasible points below the cutoff, leave the root LP as it
  // was, and (iteration limits only) return the same point on a rerun.
  int found = 0;
  for (unsigned seed = 1; seed <= 30; ++seed) {
    MipProblem p = seed % 2 ? mixed(seed) : set_cover(seed, 60, 40);
    DualSimplex lp;
    lp.load(p.A, p.cost, p.col_lb, p.col_ub, p.row_lb, p.row_ub);
    if (lp.solve(LpOptions{}) != LpStatus::Optimal) continue;
    const double obj0 = lp.objective();
    SubMipLimits lim;
    lim.node_limit = 200;
    lim.time_limit = 1e9;
    lim.iteration_limit = 200000;
    SubMipStats ss;
    lim.stats = &ss;
    std::vector<double> sol, sol2;
    LpRensStats st;
    const bool ok = lp_rens(p, lp, kInf, lim, sol, 0.7, 4, &st);
    CHECK(lp.objective() == obj0);
    for (Index j = 0; j < p.num_cols(); ++j) CHECK(lp.col_lower(j) == p.col_lb[j] && lp.col_upper(j) == p.col_ub[j]);
    if (!ok) continue;
    ++found;
    CHECK(p.max_violation(sol) <= 1e-6);
    CHECK(st.attempts >= 1);
    CHECK(ss.solves >= 1);
    // A cutoff at the found value: anything returned must be strictly better.
    const double cut = p.objective(sol);
    if (lp_rens(p, lp, cut, lim, sol2, 0.7, 4)) CHECK(p.objective(sol2) < cut);
    std::vector<double> again;
    REQUIRE(lp_rens(p, lp, kInf, lim, again, 0.7, 4));
    CHECK(again == sol);
  }
  CHECK(found >= 10);
}
