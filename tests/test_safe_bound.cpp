// Safe dual bound of the dual simplex (DualSimplex::safe_bound / safe_objective): equal to the LP
// objective on clean LPs, a valid lower bound for any row multipliers (perturbed and dual infeasible
// ones), and rounding cases where the bound computed without its error terms exceeds the optimum.
// Dual feasibility in unscaled terms (row duals of upscaled rows), and the check of infeasibility proofs
// (DualSimplex::farkas_bound): rays of infeasible LPs pass it, rays that ignore a column do not.

#include <cmath>
#include <random>
#include <vector>

#include "milp/lp/simplex.hpp"
#include "test_framework.hpp"

using namespace milp;

namespace {

struct Lp {
  SparseMatrix A;  // CSC
  std::vector<double> c, cl, cu, rl, ru;
};

// Random LP with ranged, one-sided and equality rows and boxed, one-sided and free columns.
Lp random_lp(std::mt19937& rng, Index m, Index n) {
  Lp lp;
  std::vector<Index> J, I;
  std::vector<double> V;
  for (Index j = 0; j < n; ++j)
    for (Index i = 0; i < m; ++i)
      if (rng() % 100 < 35) J.push_back(j), I.push_back(i), V.push_back(static_cast<int>(rng() % 21) - 10);
  lp.A = SparseMatrix::from_triplets(n, m, J, I, V);
  for (Index j = 0; j < n; ++j) {
    lp.c.push_back(static_cast<int>(rng() % 21) - 10);
    const int kind = static_cast<int>(rng() % 6);
    lp.cl.push_back(kind == 4 || kind == 5 ? -kInf : -static_cast<double>(rng() % 5));
    lp.cu.push_back(kind == 3 || kind == 5 ? kInf : static_cast<double>(rng() % 5));
  }
  for (Index i = 0; i < m; ++i) {
    const int kind = static_cast<int>(rng() % 5);
    const double lo = -static_cast<double>(rng() % 30), hi = static_cast<double>(rng() % 30);
    lp.rl.push_back(kind == 1 ? -kInf : kind == 3 ? hi : lo);
    lp.ru.push_back(kind == 2 ? kInf : hi);
  }
  return lp;
}

// Reference optimum: the LP's own objective, accurate to ~1e-12 on these small integer data.
double rel_tol(double v) { return 1e-7 * (1.0 + std::abs(v)); }

}  // namespace

TEST(safe_bound_equals_objective_on_clean_lps) {
  std::mt19937 rng(17);
  int optimal = 0;
  for (int t = 0; t < 400; ++t) {
    Lp lp = random_lp(rng, 1 + rng() % 25, 1 + rng() % 35);
    DualSimplex s;
    s.load(lp.A, lp.c, lp.cl, lp.cu, lp.rl, lp.ru);
    if (s.solve(LpOptions{}) != LpStatus::Optimal) continue;
    ++optimal;
    const double obj = s.objective();
    const double sb = s.safe_bound(s.row_duals(), 1e-7);
    CHECK(sb <= obj + rel_tol(obj));
    CHECK_NEAR(sb, obj, rel_tol(obj));
    bool below = true;
    CHECK_EQ(s.safe_objective(1e-7, &below), obj);
    CHECK(!below);
    // After bound changes and a warm re-solve (local column bounds enter the bound).
    for (int k = 0; k < 3; ++k) {
      const Index j = static_cast<Index>(rng() % lp.c.size());
      const double v = s.x()[j];
      if (rng() % 2) lp.cu[j] = std::max(lp.cl[j], std::floor(v));
      else lp.cl[j] = std::min(lp.cu[j], std::ceil(v));
      s.set_col_bounds(j, lp.cl[j], lp.cu[j]);
    }
    if (s.solve(LpOptions{}) != LpStatus::Optimal) continue;
    const double obj2 = s.objective();
    CHECK_NEAR(s.safe_bound(s.row_duals(), 1e-7), obj2, rel_tol(obj2));
  }
  CHECK(optimal > 100);
}

TEST(safe_bound_valid_for_perturbed_duals) {
  std::mt19937 rng(23);
  std::normal_distribution<double> gauss(0.0, 1.0);
  int dropped = 0, cases = 0;
  for (int t = 0; t < 300; ++t) {
    Lp lp = random_lp(rng, 2 + rng() % 20, 2 + rng() % 30);
    DualSimplex s;
    s.load(lp.A, lp.c, lp.cl, lp.cu, lp.rl, lp.ru);
    if (s.solve(LpOptions{}) != LpStatus::Optimal) continue;
    const double opt = s.objective();
    const std::vector<double> y0 = s.row_duals();
    for (const double scale : {1e-9, 1e-6, 1e-3, 1e-1, 1.0, 10.0}) {
      std::vector<double> y = y0;
      for (double& v : y) v += scale * gauss(rng);
      ++cases;
      const double sb = s.safe_bound(y, 1e-7);
      CHECK(sb <= opt + rel_tol(opt));
      if (scale >= 1e-1 && sb < opt - 1e-3) ++dropped;
    }
    // A dual of the wrong sign on a one-sided row (the side it needs is infinite): dropped, still valid.
    for (Index i = 0; i < static_cast<Index>(y0.size()); ++i) {
      if (std::isfinite(lp.rl[i]) && std::isfinite(lp.ru[i])) continue;
      std::vector<double> y = y0;
      y[i] = std::isfinite(lp.rl[i]) ? -0.5 : 0.5;
      ++cases;
      CHECK(s.safe_bound(y, 1e-7) <= opt + rel_tol(opt));
    }
  }
  CHECK(cases > 1000);
  CHECK(dropped > 100);  // large dual infeasibilities show as a lower bound
}

TEST(safe_bound_product_rounding) {
  // min c x s.t. x >= 3 (as a row), c = fl(1/3): the optimum 3c = 1 - 2^-54 is not a double, and
  // fl(3c) = 1.0 lies above it (tie to even). With the LP's dual y = c the bound must stay below 1.
  const double c = 1.0 / 3.0;
  DualSimplex s;
  SparseMatrix A = SparseMatrix::from_triplets(1, 1, {0}, {0}, {1.0});
  s.load(A, {c}, {0.0}, {10.0}, {3.0}, {kInf});
  REQUIRE(s.solve(LpOptions{}) == LpStatus::Optimal);
  CHECK(s.safe_bound({c}, 1e-7) < 1.0);
  CHECK(s.safe_bound(s.row_duals(), 1e-7) < 1.0);
  // The same through a column bound (r = c, x >= 3).
  DualSimplex s2;
  SparseMatrix B = SparseMatrix::from_triplets(1, 1, {0}, {0}, {1.0});
  s2.load(B, {c}, {3.0}, {10.0}, {-kInf}, {kInf});
  REQUIRE(s2.solve(LpOptions{}) == LpStatus::Optimal);
  CHECK(s2.safe_bound({0.0}, 1e-7) < 1.0);
}

TEST(safe_bound_reduced_cost_rounding) {
  // min 2^53 x1 + x2 - 2^53 x3 s.t. x0 + x1 = 0, x0 + x2 = 0, x0 + x3 = 0, x in [0, 10] x [-10, 10]^3:
  // optimum -10 (x0 = 10). With the optimal duals y = (2^53, 1, -2^53), r_0 = 0 - (2^53 + 1 - 2^53) = -1,
  // but summed in row order 2^53 + 1 rounds to 2^53 and the computed r_0 is 0: without the error term
  // of r_0 (times |u_0| = 10) the bound would come out as 0 > -10.
  const double big = std::ldexp(1.0, 53);
  SparseMatrix A = SparseMatrix::from_triplets(4, 3, {0, 0, 0, 1, 2, 3}, {0, 1, 2, 0, 1, 2}, {1, 1, 1, 1, 1, 1});
  DualSimplex s;
  s.load(A, {0.0, big, 1.0, -big}, {0.0, -10.0, -10.0, -10.0}, {10.0, 10.0, 10.0, 10.0}, {0.0, 0.0, 0.0},
         {0.0, 0.0, 0.0});
  const double sb = s.safe_bound({big, 1.0, -big}, 1e-7);
  CHECK(sb <= -10.0);
  CHECK(std::isfinite(sb));
}

TEST(safe_bound_infinite_bounds) {
  // min x0 + x1 s.t. x0 - x1 >= 2, x0 in [0, inf), x1 in [-5, inf): optimum -5 (x0 = 0, x1 = -5, y = 0).
  SparseMatrix A = SparseMatrix::from_triplets(2, 1, {0, 1}, {0, 0}, {1.0, -1.0});
  DualSimplex s;
  s.load(A, {1.0, 1.0}, {0.0, -5.0}, {kInf, kInf}, {2.0}, {kInf});
  REQUIRE(s.solve(LpOptions{}) == LpStatus::Optimal);
  CHECK_NEAR(s.objective(), -5.0, 1e-9);
  CHECK_NEAR(s.safe_bound(s.row_duals(), 1e-7), -5.0, 1e-9);
  // y = 1: r = (0, 2); x0's needed bound is l = 0, x1's l = -5: L = 2 + 0 - 10 = -8.
  CHECK_NEAR(s.safe_bound({1.0}, 1e-7), -8.0, 1e-9);
  // y = 2: r = (-1, 3): x0 needs its infinite upper bound.
  CHECK(s.safe_bound({2.0}, 1e-7) == -kInf);
  // A r_j within tau of 0 on the infinite side counts as zero (with an error term): y = 1 + 1e-9.
  const double sb = s.safe_bound({1.0 + 1e-9}, 1e-7);
  CHECK(std::isfinite(sb));
  CHECK(sb <= -5.0);
  // The wrong sign on the row (y < 0 needs the infinite upper side): dropped.
  CHECK(s.safe_bound({-1.0}, 1e-7) <= -5.0);
}

TEST(safe_bound_random_cutoff_claims) {
  // CutoffReached rests on the Lagrangian bound with error terms: a cutoff below the optimum by more
  // than the tolerances must never stop the LP, and a stop happens only for cutoffs below the optimum.
  std::mt19937 rng(29);
  for (int t = 0; t < 200; ++t) {
    Lp lp = random_lp(rng, 2 + rng() % 20, 2 + rng() % 30);
    DualSimplex s;
    s.load(lp.A, lp.c, lp.cl, lp.cu, lp.rl, lp.ru);
    if (s.solve(LpOptions{}) != LpStatus::Optimal) continue;
    const double opt = s.objective();
    for (const double delta : {-1.0, -1e-3, 1e-3, 1.0}) {
      DualSimplex f;
      f.load(lp.A, lp.c, lp.cl, lp.cu, lp.rl, lp.ru);
      LpOptions o;
      o.objective_cutoff = opt + delta;
      const LpStatus st = f.solve(o);
      if (st == LpStatus::CutoffReached) CHECK(opt > o.objective_cutoff - rel_tol(opt));
      else CHECK(st == LpStatus::Optimal);
    }
  }
}

TEST(safe_bound_catches_early_stops_on_min_max_lps) {
  // The min-max LPs of lp_negligible_entries_do_not_skew_row_scaling (test_lp.cpp): with a coefficient of
  // 8.9e-16 in a row whose largest is ~100, a row scaling that let it set the geometric scale (2^22) made the
  // dual tolerance on the row's logical accept unscaled reduced costs of ~0.1 and the LP stop above the
  // optimum. Whatever the scaling does, the safe objective must not exceed the optimum of the LP without the
  // negligible entries (the same LP up to 1e-15 x_k terms), and an early stop must show as below.
  std::mt19937 rng(97);
  int checked = 0;
  for (int t = 0; t < 60; ++t) {
    const Index rows = 5 + static_cast<Index>(rng() % 20), nx = 8 + static_cast<Index>(rng() % 30);
    const Index n = nx + 1, m = rows + 1;
    std::vector<Index> J, I, Jt, It;
    std::vector<double> V, Vt;
    for (Index i = 0; i < rows; ++i) {
      J.push_back(0), I.push_back(i), V.push_back(1.0);  // t
      for (Index k = 0; k < nx; ++k) {
        if (rng() % 100 < 60) J.push_back(1 + k), I.push_back(i), V.push_back(static_cast<double>(static_cast<int>(rng() % 4001) - 2000) / 32.0);
      }
    }
    for (Index k = 0; k < nx; ++k) J.push_back(1 + k), I.push_back(rows), V.push_back(1.0);
    Jt = J, It = I, Vt = V;
    for (Index i = 0; i < rows; ++i) {
      Index k = 1 + static_cast<Index>(rng() % nx);
      bool present = false;
      for (size_t e = 0; e < J.size(); ++e) present = present || (I[e] == i && J[e] == k);
      if (present) continue;
      Jt.push_back(k), It.push_back(i), Vt.push_back(8.881784197001252e-16);
    }
    std::vector<double> c(n, 0.0), cl(n, 0.0), cu(n, 1.0), rl(m, 0.0), ru(m, kInf);
    c[0] = 1.0, cl[0] = -kInf, cu[0] = kInf;
    rl[rows] = ru[rows] = static_cast<double>(1 + rng() % (nx / 2));
    DualSimplex a, b;
    a.load(SparseMatrix::from_triplets(n, m, J, I, V), c, cl, cu, rl, ru);
    b.load(SparseMatrix::from_triplets(n, m, Jt, It, Vt), c, cl, cu, rl, ru);
    REQUIRE(a.solve(LpOptions{}) == LpStatus::Optimal);
    REQUIRE(b.solve(LpOptions{}) == LpStatus::Optimal);
    const double opt = a.objective(), slack = static_cast<double>(nx) * 1e-15;
    bool below = false;
    const double safe = b.safe_objective(1e-7, &below);
    CHECK(safe <= opt + slack + rel_tol(opt));
    if (b.objective() > opt + rel_tol(opt)) CHECK(below);
    ++checked;
  }
  CHECK(checked == 60);
}

TEST(lp_dual_feasible_in_unscaled_terms) {
  // min 1e-3 x0 s.t. 100 x0 + 1e-6 x1 <= 50, x in [0, 1]^2: optimum 0 (x0 = 0). The row's scale is 2^7 (geometric
  // mean 0.01): 12800 x0 + 1.28e-4 x1 <= 6400. From the basis {x0} with the row tight (x0 = 0.5, objective 5e-4)
  // the row's scaled dual is 1e-3 / 12800 = 7.8e-8 of the wrong sign, inside the 1e-7 tolerance in scaled terms
  // but 1e-5 unscaled: the LP must go on to the optimum, and its safe bound then equals the objective.
  SparseMatrix A = SparseMatrix::from_triplets(2, 1, {0, 1}, {0, 0}, {100.0, 1e-6});
  DualSimplex s;
  s.load(A, {1e-3, 0.0}, {0.0, 0.0}, {1.0, 1.0}, {-kInf}, {50.0});
  s.set_basis(Basis{{VarStatus::Basic, VarStatus::AtLower, VarStatus::AtLower}});
  REQUIRE(s.solve(LpOptions{}) == LpStatus::Optimal);
  CHECK(std::abs(s.objective()) <= 1e-9);
  bool below = true;
  CHECK_EQ(s.safe_objective(1e-7, &below), s.objective());
  CHECK(!below);
}

TEST(lp_dual_feasible_in_unscaled_terms_min_max) {
  // The min-max LPs of safe_bound_catches_early_stops_on_min_max_lps with a small entry of 1e-6 (above the row
  // scaling's 1e-9 cut) in every min-max row: row scales 2^6 - 2^7, scaled entries up to ~8000. At optimality
  // every nonbasic variable's unscaled reduced cost (logicals: minus the row's dual) has the right sign within
  // the 1e-7 tolerance, the safe bound equals the objective, and the objective is that of the LP without the
  // small entries (up to their effect, nx 1e-6).
  std::mt19937 rng(41);
  int checked = 0;
  for (int t = 0; t < 60; ++t) {
    const Index rows = 5 + static_cast<Index>(rng() % 20), nx = 8 + static_cast<Index>(rng() % 30);
    const Index n = nx + 1, m = rows + 1;
    std::vector<Index> J, I, Jt, It;
    std::vector<double> V, Vt;
    for (Index i = 0; i < rows; ++i) {
      J.push_back(0), I.push_back(i), V.push_back(1.0);  // t
      for (Index k = 0; k < nx; ++k) {
        if (rng() % 100 < 60) J.push_back(1 + k), I.push_back(i), V.push_back(static_cast<double>(static_cast<int>(rng() % 4001) - 2000) / 32.0);
      }
    }
    for (Index k = 0; k < nx; ++k) J.push_back(1 + k), I.push_back(rows), V.push_back(1.0);
    Jt = J, It = I, Vt = V;
    for (Index i = 0; i < rows; ++i) {
      Index k = 1 + static_cast<Index>(rng() % nx);
      bool present = false;
      for (size_t e = 0; e < J.size(); ++e) present = present || (I[e] == i && J[e] == k);
      if (present) continue;
      Jt.push_back(k), It.push_back(i), Vt.push_back(1e-6);
    }
    std::vector<double> c(n, 0.0), cl(n, 0.0), cu(n, 1.0), rl(m, 0.0), ru(m, kInf);
    c[0] = 1.0, cl[0] = -kInf, cu[0] = kInf;
    rl[rows] = ru[rows] = static_cast<double>(1 + rng() % (nx / 2));
    DualSimplex a, b;
    a.load(SparseMatrix::from_triplets(n, m, J, I, V), c, cl, cu, rl, ru);
    b.load(SparseMatrix::from_triplets(n, m, Jt, It, Vt), c, cl, cu, rl, ru);
    REQUIRE(a.solve(LpOptions{}) == LpStatus::Optimal);
    REQUIRE(b.solve(LpOptions{}) == LpStatus::Optimal);
    bool ok = true;
    for (Index j = 0; j < n + m; ++j) {
      const VarStatus st = b.status(j);
      if (st == VarStatus::Basic || b.var_lower(j) == b.var_upper(j)) continue;
      const double d = b.reduced_cost(j);
      if (st == VarStatus::AtLower) ok = ok && d >= -1e-7;
      else if (st == VarStatus::AtUpper) ok = ok && d <= 1e-7;
      else ok = ok && std::abs(d) <= 1e-7;
    }
    CHECK(ok);
    bool below = true;
    CHECK_EQ(b.safe_objective(1e-7, &below), b.objective());
    CHECK(!below);
    const double slack = static_cast<double>(nx) * 1e-6;
    if (!CHECK(std::abs(b.objective() - a.objective()) <= slack + rel_tol(a.objective()))) {
      std::printf("  case %d: %.10g without, %.10g with the small entries\n", t, a.objective(), b.objective());
    }
    ++checked;
  }
  CHECK(checked == 60);
}

namespace {
// Zero-cost Lagrangian of the last solve's ray, with the sign that proves infeasibility (farkas_proof does not
// normalize it) and DualSimplex's tau (1e-9 of the largest multiplier).
double ray_bound(const DualSimplex& s) {
  std::vector<double> y;
  if (!s.farkas_proof(y)) return -kInf;
  double big = 0.0;
  for (double v : y) big = std::max(big, std::abs(v));
  std::vector<double> ny = y;
  for (double& v : ny) v = -v;
  return std::max(s.farkas_bound(y, 1e-9 * big), s.farkas_bound(ny, 1e-9 * big));
}
}  // namespace

TEST(farkas_proofs_verified_on_infeasible_lps) {
  // Random LPs, rows multiplied by 1e-4 or 1e3 (row scales other than 1), solved cold and warm after column
  // bound changes (an updated factorization): every Infeasible verdict comes with a ray whose zero-cost
  // Lagrangian bound on the unscaled data is positive, and no ray fails the solver's own check.
  std::mt19937 rng(53);
  const DualSimplex::CheckCounts c0 = DualSimplex::check_counts();
  int cold = 0, warm = 0;
  for (int t = 0; t < 400; ++t) {
    Lp lp = random_lp(rng, 2 + rng() % 25, 2 + rng() % 35);
    std::vector<double> f(lp.rl.size(), 1.0);
    for (double& v : f) {
      const int k = static_cast<int>(rng() % 3);
      v = k == 0 ? 1e-4 : k == 1 ? 1e3 : 1.0;
    }
    for (Index j = 0; j < static_cast<Index>(lp.c.size()); ++j)
      for (Index p = lp.A.begin(j); p < lp.A.end(j); ++p) lp.A.value[p] *= f[lp.A.index[p]];
    for (size_t i = 0; i < f.size(); ++i) lp.rl[i] *= f[i], lp.ru[i] *= f[i];
    DualSimplex s;
    s.load(lp.A, lp.c, lp.cl, lp.cu, lp.rl, lp.ru);
    LpStatus st = s.solve(LpOptions{});
    if (st == LpStatus::Infeasible) {
      ++cold;
      CHECK(ray_bound(s) > 0.0);
      continue;
    }
    // Shrink column boxes away from the LP point until the LP becomes infeasible.
    for (int k = 0; k < 12 && st == LpStatus::Optimal; ++k) {
      const Index j = static_cast<Index>(rng() % lp.c.size());
      const double v = s.x()[j];
      if (rng() % 2) lp.cu[j] = std::max(lp.cl[j], std::floor(v) - 1.0);
      else lp.cl[j] = std::min(lp.cu[j], std::ceil(v) + 1.0);
      s.set_col_bounds(j, lp.cl[j], lp.cu[j]);
      st = s.solve(LpOptions{});
      if (st == LpStatus::Infeasible) {
        ++warm;
        CHECK(ray_bound(s) > 0.0);
      }
    }
  }
  CHECK(cold > 20);
  CHECK(warm > 20);
  const DualSimplex::CheckCounts c1 = DualSimplex::check_counts();
  if (DualSimplex::farkas_check_enabled()) CHECK(c1.farkas_checked - c0.farkas_checked >= cold + warm);
  CHECK_EQ(c1.farkas_failed - c0.farkas_failed, 0);
}

TEST(farkas_check_rejects_ignored_columns) {
  // x + 5e-10 z >= 1, x in [0, 0.5], z in [0, 1e10]: feasible (z = 1e9), but the dual ratio test ignores z's
  // entry (below 1e-9; nor does it count in the row scale) and finds a ray that proves infeasibility only
  // without z: y = 1 gives 1 - 0.5 - 5e-10 * 1e10 < 0. The solve must not report Infeasible.
  SparseMatrix A = SparseMatrix::from_triplets(2, 1, {0, 1}, {0, 0}, {1.0, 5e-10});
  DualSimplex s;
  s.load(A, {1.0, 0.0}, {0.0, 0.0}, {0.5, 1e10}, {1.0}, {kInf});
  CHECK(s.solve(LpOptions{}) != LpStatus::Infeasible);
  CHECK(s.farkas_bound({1.0}, 1e-9) <= 0.0);
  // With z <= 1e8 the LP is infeasible (x + 5e-10 z <= 0.55) and the same ray proves it.
  s.set_col_bounds(1, 0.0, 1e8);
  s.reset_basis();
  REQUIRE(s.solve(LpOptions{}) == LpStatus::Infeasible);
  CHECK(ray_bound(s) > 0.0);
  CHECK_NEAR(s.farkas_bound({1.0}, 1e-9), 0.45, 1e-9);
}

TEST(farkas_bound_valid_on_feasible_lps) {
  // No multipliers prove a feasible LP infeasible: the zero-cost Lagrangian stays <= 0 for random y.
  std::mt19937 rng(59);
  std::normal_distribution<double> gauss(0.0, 1.0);
  int cases = 0;
  for (int t = 0; t < 200; ++t) {
    Lp lp = random_lp(rng, 2 + rng() % 20, 2 + rng() % 30);
    DualSimplex s;
    s.load(lp.A, lp.c, lp.cl, lp.cu, lp.rl, lp.ru);
    if (s.solve(LpOptions{}) != LpStatus::Optimal) continue;
    for (int k = 0; k < 10; ++k) {
      std::vector<double> y(lp.rl.size());
      for (double& v : y) v = rng() % 3 == 0 ? 0.0 : gauss(rng);
      ++cases;
      CHECK(s.farkas_bound(y, 1e-9) <= 0.0);
    }
  }
  CHECK(cases > 500);
}
