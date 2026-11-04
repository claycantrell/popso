// LP duals and fixed-MIP duals (solver/duals.hpp, SolverParams::duals): KKT conditions on random
// LPs (presolve on and off, minimization and maximization) and on the generated LP instances; the
// fixed-MIP duals against an LP solve with the integers fixed.

#include <algorithm>
#include <atomic>
#include <cmath>
#include <random>
#include <string>
#include <vector>

#include "milp/core/interrupt.hpp"
#include "milp/core/model.hpp"
#include "milp/io/mps_reader.hpp"
#include "milp/solver/duals.hpp"
#include "milp/solver/solver.hpp"
#include "test_framework.hpp"

using namespace milp;

namespace {

SolverParams quiet_params(int duals) {
  SolverParams sp;
  sp.verbosity = 0;
  sp.threads = 1;
  sp.duals = duals;
  return sp;
}

// Scale of the data for the relative tolerances.
double data_scale(const Model& m) {
  double s = 1.0;
  for (double c : m.obj) s = std::max(s, std::abs(c));
  for (double v : m.A.value) s = std::max(s, std::abs(v));
  return s;
}

struct Kkt {
  double stationarity = 0.0;    // |d - (c - A^T y)|
  double sign = 0.0;            // multipliers pushing against an infinite bound
  double complementarity = 0.0;  // multiplier x distance to its bound
  double max() const { return std::max({stationarity, sign, complementarity}); }
};

// Multiplier z of a value v in [lo, hi], minimization form (z >= 0 pushes against lo, z <= 0
// against hi): wrong sign where the bound is infinite, and z times the distance to the bound.
void check_side(double v, double lo, double hi, double z, double scale, Kkt& k) {
  const double zp = std::max(z, 0.0), zm = std::max(-z, 0.0);
  if (lo == hi) return;  // fixed: any sign
  if (!is_finite(lo)) k.sign = std::max(k.sign, zp);
  else k.complementarity = std::max(k.complementarity, zp * std::abs(v - lo) / scale);
  if (!is_finite(hi)) k.sign = std::max(k.sign, zm);
  else k.complementarity = std::max(k.complementarity, zm * std::abs(hi - v) / scale);
}

// KKT conditions of an LP result with (x, y, d) in the model's sense: minimization y_i >= 0 at
// an active row lower bound, d_j >= 0 at a column lower bound; maximization the opposite.
// col_lb / col_ub override the model's bounds (fixed-MIP duals).
Kkt lp_kkt(const Model& m, const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& d,
           const std::vector<double>& col_lb, const std::vector<double>& col_ub) {
  const double sign = m.sense == ObjSense::Maximize ? -1.0 : 1.0;
  const double scale = data_scale(m);
  Kkt k;
  std::vector<double> act(m.num_rows(), 0.0);
  for (Index j = 0; j < m.num_cols(); ++j) {
    double g = m.obj[j];
    for (Index p = m.A.begin(j); p < m.A.end(j); ++p) {
      g -= m.A.value[p] * y[m.A.index[p]];
      act[m.A.index[p]] += m.A.value[p] * x[j];
    }
    k.stationarity = std::max(k.stationarity, std::abs(g - d[j]) / scale);
    check_side(x[j], col_lb[j], col_ub[j], sign * d[j], 1.0 + std::abs(x[j]), k);
  }
  for (Index i = 0; i < m.num_rows(); ++i) check_side(act[i], m.row_lb[i], m.row_ub[i], sign * y[i], 1.0 + std::abs(act[i]), k);
  return k;
}

Kkt lp_kkt(const Model& m, const SolveResult& r) { return lp_kkt(m, r.x, r.row_dual, r.reduced_cost, m.col_lb, m.col_ub); }

// Dual objective of (y, d) for the model with these column bounds (model's sense, incl. offset):
// sum over rows and columns of the multiplier times the bound it pushes against.
double dual_objective(const Model& m, const std::vector<double>& y, const std::vector<double>& d,
                      const std::vector<double>& col_lb, const std::vector<double>& col_ub) {
  const double sign = m.sense == ObjSense::Maximize ? -1.0 : 1.0;
  double v = 0.0;
  auto term = [&](double z, double lo, double hi) {
    const double zs = sign * z;  // minimization form
    if (zs > 0.0) return is_finite(lo) ? zs * lo : (zs > 1e-9 ? -kInf : 0.0);
    if (zs < 0.0) return is_finite(hi) ? zs * hi : (zs < -1e-9 ? -kInf : 0.0);
    return 0.0;
  };
  for (Index i = 0; i < m.num_rows(); ++i) v += term(y[i], m.row_lb[i], m.row_ub[i]);
  for (Index j = 0; j < m.num_cols(); ++j) v += term(d[j], col_lb[j], col_ub[j]);
  return sign * v + m.obj_offset;
}

// A random LP around a feasible point x0 with every kind of column and row bound, plus
// structures the presolve reduces (singleton rows, duplicate and parallel rows and columns,
// doubleton equations, fixed and empty columns, free rows).
Model random_lp(std::mt19937& rng, Index n, Index m, bool integers) {
  std::uniform_real_distribution<double> u(-1.0, 1.0);
  ModelBuilder b;
  std::vector<double> x0(n);
  std::vector<double> lo(n), hi(n);
  for (Index j = 0; j < n; ++j) {
    x0[j] = std::round(4.0 * u(rng));
    const int kind = static_cast<int>(rng() % 7);
    lo[j] = x0[j] - (rng() % 3);
    hi[j] = x0[j] + (rng() % 3);
    if (kind == 0) lo[j] = -kInf;
    if (kind == 1) hi[j] = kInf;
    if (kind == 2) lo[j] = hi[j] = x0[j];  // fixed
    if (kind == 3) lo[j] = -kInf, hi[j] = kInf;
    const bool is_int = integers && rng() % 2 == 0;
    // Free columns get a cost of either sign only if the rows bound them; keep the LP bounded by
    // boxing them in a far row pair below.
    b.add_col("x" + std::to_string(j), std::round(10.0 * u(rng)) / 2.0, lo[j], hi[j],
              is_int ? VarType::Integer : VarType::Continuous);
  }
  // Duplicate one column (parallel columns with proportional costs).
  if (n >= 2 && rng() % 2 == 0) {
    Model& md = b.model();
    md.obj[1] = 2.0 * md.obj[0];
  }
  std::vector<std::vector<double>> rows;
  for (Index i = 0; i < m; ++i) {
    std::vector<double> a(n, 0.0);
    const int kind = static_cast<int>(rng() % 8);
    if (kind == 0) {
      a[rng() % n] = std::round(4.0 * u(rng)) + 0.5;  // singleton row
    } else if (kind == 1 && n >= 2) {
      const Index j1 = static_cast<Index>(rng() % n), j2 = (j1 + 1) % n;  // doubleton
      a[j1] = 1.0;
      a[j2] = std::round(3.0 * u(rng)) + 0.25;
    } else if (kind == 2 && !rows.empty()) {
      a = rows[rng() % rows.size()];  // parallel row
      const double f = rng() % 2 ? 2.0 : -0.5;
      for (double& v : a) v *= f;
    } else {
      for (Index j = 0; j < n; ++j) a[j] = rng() % 3 == 0 ? 0.0 : std::round(8.0 * u(rng)) / 4.0;
    }
    if (n >= 2 && rng() % 4 == 0) a[1] = 2.0 * a[0];  // keeps column 1 parallel to column 0 sometimes
    double act = 0.0;
    for (Index j = 0; j < n; ++j) act += a[j] * x0[j];
    double rl = act - (rng() % 3), rh = act + (rng() % 3);
    const int side = static_cast<int>(rng() % 6);
    if (side == 0) rl = -kInf;
    if (side == 1) rh = kInf;
    if (side == 2) rl = rh = act;  // equality
    if (side == 3 && rng() % 4 == 0) rl = -kInf, rh = kInf;  // free row
    const Index r = b.add_row("r" + std::to_string(i), rl, rh);
    for (Index j = 0; j < n; ++j)
      if (a[j] != 0.0) b.add_coef(r, j, a[j]);
    rows.push_back(a);
  }
  // Box rows keep the LP bounded: -B <= x_j <= B as rows (not bounds) for the unbounded columns.
  for (Index j = 0; j < n; ++j) {
    if (is_finite(lo[j]) && is_finite(hi[j])) continue;
    const Index r = b.add_row("box" + std::to_string(j), -20.0, 20.0);
    b.add_coef(r, j, 1.0);
    if (n >= 2 && rng() % 2 == 0) b.add_coef(r, (j + 1) % n, 0.5);
  }
  Model md = b.build();
  md.obj_offset = std::round(4.0 * u(rng));
  return md;
}

Model maximized(const Model& md) {
  Model mx = md;
  mx.sense = ObjSense::Maximize;
  for (double& c : mx.obj) c = -c;
  mx.obj_offset = -mx.obj_offset;
  return mx;
}

// Basic statuses: one basic variable per row.
bool basis_square(const SolveResult& r) {
  long long basic = 0;
  for (int8_t s : r.col_basis) basic += s == 0;
  for (int8_t s : r.row_basis) basic += s == 0;
  return basic == static_cast<long long>(r.row_basis.size());
}

}  // namespace

TEST(lp_duals_kkt_random) {
  std::mt19937 rng(91);
  int checked = 0, presolved = 0;
  for (int trial = 0; trial < 150; ++trial) {
    const Index n = 2 + static_cast<Index>(rng() % 9), m = 1 + static_cast<Index>(rng() % 8);
    const Model md = random_lp(rng, n, m, false);
    for (int sense = 0; sense < 2; ++sense) {
      const Model mm = sense ? maximized(md) : md;
      double ref = kInf;
      for (int presolve = 1; presolve >= 0; --presolve) {
        SolverParams sp = quiet_params(1);
        sp.presolve = presolve != 0;
        const SolveResult r = solve(mm, sp);
        if (r.status != Status::Optimal) {
          CHECK(r.status == Status::Infeasible || r.status == Status::InfeasibleOrUnbounded);
          CHECK(!r.has_duals);
          continue;
        }
        if (!is_finite(ref)) ref = r.objective;
        CHECK_NEAR(r.objective, ref, 1e-7 * (1.0 + std::abs(ref)));
        REQUIRE(r.has_duals);
        REQUIRE(r.row_dual.size() == static_cast<size_t>(mm.num_rows()));
        REQUIRE(r.reduced_cost.size() == static_cast<size_t>(mm.num_cols()));
        CHECK(!r.fixed_mip_duals);
        const Kkt k = lp_kkt(mm, r);
        CHECK(k.max() <= 1e-7);
        if (k.max() > 1e-7) {
          std::printf("  trial %d sense %d presolve %d: stationarity %.3g sign %.3g compl %.3g\n", trial, sense, presolve,
                      k.stationarity, k.sign, k.complementarity);
        }
        // Strong duality.
        const double dobj = dual_objective(mm, r.row_dual, r.reduced_cost, mm.col_lb, mm.col_ub);
        CHECK_NEAR(dobj, r.objective, 1e-7 * (1.0 + std::abs(r.objective)));
        CHECK(basis_square(r));
        ++checked;
        presolved += presolve;
      }
    }
  }
  CHECK(checked >= 300);
  CHECK(presolved >= 150);
}

TEST(lp_duals_off_and_qp_unchanged) {
  std::mt19937 rng(5);
  const Model md = random_lp(rng, 5, 4, false);
  const SolveResult r = solve(md, quiet_params(0));
  CHECK(!r.has_duals && r.row_dual.empty() && r.col_basis.empty());
  // A MILP gets fixed-MIP duals only with duals = 2.
  Model mi = random_lp(rng, 6, 4, true);
  const SolveResult ri = solve(mi, quiet_params(1));
  if (mi.num_integers() > 0) CHECK(!ri.has_duals);
}

TEST(lp_duals_generated_instances) {
  const std::string dir = std::string(MILP_SOURCE_DIR) + "/benchmarks/instances/";
  const char* names[] = {"randlp_small_1", "randlp_medium_1", "transport_small_1", "transport_medium_1",
                         "lpsanity_beale_cycling_0", "lpsanity_degenerate_0", "lpsanity_dual_degenerate_0",
                         "lpsanity_empty_col_0", "lpsanity_empty_row_0", "lpsanity_extra_free_row_0",
                         "lpsanity_feasible_0", "lpsanity_fixed_0", "lpsanity_fixed_format_0", "lpsanity_freevars_0",
                         "lpsanity_illscaled_0", "lpsanity_klee_minty_0", "lpsanity_negative_lb_0",
                         "lpsanity_negative_ub_0", "lpsanity_objoffset_0", "lpsanity_objsense_inline_0",
                         "lpsanity_objsense_maximize_0", "lpsanity_ranges_0", "lpsanity_redundant_eq_0",
                         "lpsanity_upper_only_0"};
  int checked = 0;
  for (const char* name : names) {
    Model md;
    try {
      md = read_mps(dir + name + ".mps");
    } catch (const std::exception&) {
      continue;  // instance files are optional (generated)
    }
    REQUIRE(md.num_integers() == 0);
    for (int presolve = 0; presolve < 2; ++presolve) {
      SolverParams sp = quiet_params(1);
      sp.threads = 2;
      sp.presolve = presolve != 0;
      const SolveResult r = solve(md, sp);
      CHECK(r.status == Status::Optimal);
      if (r.status != Status::Optimal) continue;
      CHECK(r.has_duals);
      if (!r.has_duals) continue;
      const Kkt k = lp_kkt(md, r);
      // Larger instances: the simplex tolerances (1e-7 primal / dual, relative to the scaling).
      CHECK(k.max() <= 1e-6);
      if (k.max() > 1e-6) {
        std::printf("  %s presolve %d: stationarity %.3g sign %.3g compl %.3g\n", name, presolve, k.stationarity, k.sign,
                    k.complementarity);
      }
      const double dobj = dual_objective(md, r.row_dual, r.reduced_cost, md.col_lb, md.col_ub);
      CHECK_NEAR(dobj, r.objective, 1e-6 * (1.0 + std::abs(r.objective)));
      CHECK(basis_square(r));
      ++checked;
    }
  }
  CHECK(checked >= 40);
}

TEST(lp_duals_long_objective_row) {
  // An LP whose objective is defined by a long equation z = sum_k c_k x_k (1200 entries, z a free column
  // singleton with cost 1, as blp-ic98 / leo1): presolve substitutes z out of the row; the duals are
  // recovered in the original space and must satisfy the KKT conditions with presolve on and off.
  std::mt19937 rng(2611);
  for (int trial = 0; trial < 4; ++trial) {
    ModelBuilder b;
    const int n = 1200;
    std::vector<Index> x;
    for (int k = 0; k < n; ++k) {
      x.push_back(b.add_col("x" + std::to_string(k), 0.0, 0.0, 1.0 + (rng() % 3), VarType::Continuous));
    }
    const Index z = b.add_col("z", 1.0, -kInf, kInf, VarType::Continuous);
    const Index obj = b.add_row("obj", 0.0, 0.0);
    for (int k = 0; k < n; ++k) {
      const double c = 1.0 + static_cast<double>(rng() % 9) - (trial % 2 ? 4.0 : 0.0);
      b.add_coef(obj, x[k], c == 0.0 ? 0.5 : c);
    }
    b.add_coef(obj, z, -1.0);
    for (int i = 0; i < 30; ++i) {  // covering rows
      const Index r = b.add_row("c" + std::to_string(i), 1.0 + (rng() % 4), kInf);
      for (int t = 0; t < 40; ++t) b.add_coef(r, x[(i * 37 + t * 29) % n], 1.0 + (rng() % 3));
    }
    Model md = b.build();
    double ref = kInf;
    for (int presolve = 1; presolve >= 0; --presolve) {
      SolverParams sp = quiet_params(1);
      sp.presolve = presolve != 0;
      const SolveResult r = solve(md, sp);
      REQUIRE(r.status == Status::Optimal);
      if (!is_finite(ref)) ref = r.objective;
      CHECK_NEAR(r.objective, ref, 1e-7 * (1.0 + std::abs(ref)));
      REQUIRE(r.has_duals);
      const Kkt k = lp_kkt(md, r);
      CHECK(k.max() <= 1e-7);
      if (k.max() > 1e-7) {
        std::printf("  trial %d presolve %d: stationarity %.3g sign %.3g compl %.3g\n", trial, presolve,
                    k.stationarity, k.sign, k.complementarity);
      }
      CHECK_NEAR(dual_objective(md, r.row_dual, r.reduced_cost, md.col_lb, md.col_ub), r.objective,
                 1e-7 * (1.0 + std::abs(r.objective)));
      CHECK(basis_square(r));
    }
  }
}

TEST(fixed_mip_duals_match_fixed_lp) {
  std::mt19937 rng(17);
  int checked = 0;
  for (int trial = 0; trial < 80; ++trial) {
    const Index n = 3 + static_cast<Index>(rng() % 8), m = 2 + static_cast<Index>(rng() % 6);
    const Model base = random_lp(rng, n, m, true);
    if (base.num_integers() == 0) continue;
    const Model md = trial % 2 ? maximized(base) : base;
    SolverParams sp = quiet_params(2);
    sp.mip_rel_gap = 0.0;
    sp.presolve = trial % 3 != 0;
    const SolveResult r = solve(md, sp);
    if (!r.has_solution) {
      CHECK(!r.has_duals);
      continue;
    }
    REQUIRE(r.has_duals);
    CHECK(r.fixed_mip_duals);
    // The LP with the integers fixed at the solution.
    Model fx = md;
    for (Index j = 0; j < md.num_cols(); ++j) {
      if (md.col_type[j] == VarType::Continuous) continue;
      fx.col_type[j] = VarType::Continuous;
      fx.col_lb[j] = fx.col_ub[j] = std::round(r.x[j]);
    }
    const SolveResult rl = solve(fx, quiet_params(1));
    REQUIRE(rl.status == Status::Optimal);
    // Stationarity, signs and strong duality of the fixed-MIP duals on the fixed LP, against
    // the fixed LP's own optimum and solution.
    const Kkt k = lp_kkt(fx, rl.x, r.row_dual, r.reduced_cost, fx.col_lb, fx.col_ub);
    CHECK(k.max() <= 1e-7);
    if (k.max() > 1e-7) {
      std::printf("  trial %d: stationarity %.3g sign %.3g compl %.3g\n", trial, k.stationarity, k.sign,
                  k.complementarity);
    }
    const double dobj = dual_objective(fx, r.row_dual, r.reduced_cost, fx.col_lb, fx.col_ub);
    CHECK_NEAR(dobj, rl.objective, 1e-7 * (1.0 + std::abs(rl.objective)));
    // Solved to optimality: the incumbent's continuous part is optimal for the fixed LP.
    if (r.status == Status::Optimal) CHECK_NEAR(rl.objective, r.objective, 1e-6 * (1.0 + std::abs(r.objective)));
    // The fixed LP's own duals have the same dual objective.
    REQUIRE(rl.has_duals);
    CHECK_NEAR(dual_objective(fx, rl.row_dual, rl.reduced_cost, fx.col_lb, fx.col_ub), rl.objective,
               1e-7 * (1.0 + std::abs(rl.objective)));
    ++checked;
  }
  CHECK(checked >= 30);
}

TEST(dual_objective_check) {
  // LP: the basis must reproduce the optimum; fixed-MIP: it may be better than the incumbent.
  CHECK(check_dual_objective(10.0, 10.0 + 1e-9, false) == DualObjectiveCheck::Consistent);
  CHECK(check_dual_objective(1e6, 1e6 * (1.0 + 1e-8), false) == DualObjectiveCheck::Consistent);
  CHECK(check_dual_objective(9.0, 10.0, false) == DualObjectiveCheck::Mismatch);
  CHECK(check_dual_objective(11.0, 10.0, false) == DualObjectiveCheck::Mismatch);
  CHECK(check_dual_objective(9.0, 10.0, true) == DualObjectiveCheck::BetterThanSolution);
  CHECK(check_dual_objective(11.0, 10.0, true) == DualObjectiveCheck::Mismatch);
  CHECK(check_dual_objective(std::nan(""), 10.0, false) == DualObjectiveCheck::Mismatch);
}

TEST(dual_recovery_is_interruptible) {
  // An interrupt or an expired time limit gives up on the duals (the crossover and the simplex
  // cleanup both stop); otherwise the same call recovers them.
  std::mt19937 rng(5);
  int checked = 0;
  for (int trial = 0; trial < 20 && checked < 5; ++trial) {
    const Model md = random_lp(rng, 8, 5, false);
    const SolveResult r = solve(md, quiet_params(1));
    if (r.status != Status::Optimal) continue;
    MipProblem p;
    REQUIRE(MipProblem::from_model(md, &p));
    CHECK(recover_lp_duals(p, {}, {}, r.x, 60.0).ok);
    CHECK(!recover_lp_duals(p, {}, {}, r.x, 0.0).ok);
    std::atomic<bool> flag{true};
    {
      InterruptScope scope(&flag);
      CHECK(!recover_lp_duals(p, {}, {}, r.x, 60.0).ok);
    }
    ++checked;
  }
  CHECK(checked == 5);
}
