// The second-order cone interior-point method (conic/socp.hpp): known optima, certificates of infeasibility and
// unboundedness.

#include <cmath>
#include <vector>

#include "milp/conic/socp.hpp"
#include "test_framework.hpp"

using namespace milp;

namespace {

// Builds a ConicProblem from dense rows (A: p x n, G: m x n).
ConicProblem make(const std::vector<double>& c, const std::vector<std::vector<double>>& A, const std::vector<double>& b,
                  const std::vector<std::vector<double>>& G, const std::vector<double>& h, Index l,
                  const std::vector<Index>& cones) {
  ConicProblem p;
  p.n = static_cast<Index>(c.size());
  p.c = c;
  p.b = b;
  p.h = h;
  p.l = l;
  p.cones = cones;
  auto csc = [&](const std::vector<std::vector<double>>& M) {
    std::vector<Index> J, I;
    std::vector<double> V;
    for (size_t i = 0; i < M.size(); ++i) {
      for (size_t j = 0; j < M[i].size(); ++j) {
        if (M[i][j] != 0.0) J.push_back(static_cast<Index>(j)), I.push_back(static_cast<Index>(i)), V.push_back(M[i][j]);
      }
    }
    return SparseMatrix::from_triplets(p.n, static_cast<Index>(M.size()), J, I, V);
  };
  p.A = csc(A);
  p.G = csc(G);
  return p;
}

}  // namespace

TEST(conic_lp_and_cone_optima) {
  ConicOptions opt;
  // LP: min -x - 2y s.t. x + y <= 1, x, y >= 0 (as G rows): -2 at (0, 1).
  {
    const ConicProblem p = make({-1, -2}, {}, {}, {{1, 1}, {-1, 0}, {0, -1}}, {1, 0, 0}, 3, {});
    const ConicResult r = solve_conic(p, opt);
    CHECK(r.status == ConicStatus::Optimal);
    CHECK_NEAR(r.pobj, -2.0, 1e-7);
    CHECK_NEAR(r.dobj, -2.0, 1e-7);
    CHECK_NEAR(r.x[1], 1.0, 1e-6);
  }
  // SOC: min t s.t. ||(x, y)|| <= t, x + y = 2 (A): sqrt(2) at x = y = 1. Variables (t, x, y); s = (t, x, y) = -G v.
  {
    const ConicProblem p = make({1, 0, 0}, {{0, 1, 1}}, {2}, {{-1, 0, 0}, {0, -1, 0}, {0, 0, -1}}, {0, 0, 0}, 0, {3});
    const ConicResult r = solve_conic(p, opt);
    CHECK(r.status == ConicStatus::Optimal);
    CHECK_NEAR(r.pobj, std::sqrt(2.0), 1e-7);
    CHECK_NEAR(r.x[1], 1.0, 1e-6);
    CHECK_NEAR(r.x[2], 1.0, 1e-6);
  }
  // Two cones and the orthant: min x + y s.t. ||(1, x - 2)|| <= y (cone 1), ||(x, y)|| <= 3 (cone 2), x >= 0.
  // On cone 1's boundary y = sqrt(1 + (x - 2)^2); x + y minimal where 1 + (x - 2) / y = 0 -> x - 2 = -y ->
  // y^2 = 1 + y^2 impossible: the minimum is at x = 0 (the orthant), y = sqrt(5).
  {
    // variables (x, y); cone 1 rows: s = (y, 1, x - 2) -> G = (0,-1), (0,0), (-1,0); h = (0, 1, -2).
    // cone 2 rows: s = (3, x, y) -> G = (0,0), (-1,0), (0,-1); h = (3, 0, 0). Orthant: s = x -> G = (-1, 0), h = 0.
    const ConicProblem p = make({1, 1}, {}, {},
                                {{-1, 0}, {0, -1}, {0, 0}, {-1, 0}, {0, 0}, {-1, 0}, {0, -1}},
                                {0, 0, 1, -2, 3, 0, 0}, 1, {3, 3});
    const ConicResult r = solve_conic(p, opt);
    CHECK(r.status == ConicStatus::Optimal);
    CHECK_NEAR(r.pobj, std::sqrt(5.0), 1e-7);
    CHECK_NEAR(r.x[0], 0.0, 1e-6);
  }
}

TEST(conic_certificates) {
  ConicOptions opt;
  // Infeasible: x >= 1 and x <= 0.
  {
    const ConicProblem p = make({1}, {}, {}, {{-1}, {1}}, {-1, 0}, 2, {});
    CHECK(solve_conic(p, opt).status == ConicStatus::PrimalInfeasible);
  }
  // Infeasible cone: ||(x, 1)|| <= t, t <= 0.5 (variables t, x).
  {
    const ConicProblem p = make({0, 0}, {}, {}, {{1, 0}, {-1, 0}, {0, 0}, {0, -1}}, {0.5, 0, 1, 0}, 1, {3});
    CHECK(solve_conic(p, opt).status == ConicStatus::PrimalInfeasible);
  }
  // Unbounded: min -x s.t. ||y|| <= x (variables x, y).
  {
    const ConicProblem p = make({-1, 0}, {}, {}, {{-1, 0}, {0, -1}}, {0, 0}, 0, {2});
    CHECK(solve_conic(p, opt).status == ConicStatus::DualInfeasible);
  }
}

TEST(conic_quadratic_objective) {
  ConicOptions opt;
  // min 1/2 (x^2 + y^2) - x - 2y s.t. x + y <= 1: on x + y = 1 the objective is x^2 - 1.5: -1.5 at (0, 1).
  ConicProblem p = make({-1, -2}, {}, {}, {{1, 1}}, {1}, 1, {});
  p.P = SparseMatrix::from_triplets(2, 2, {0, 1}, {0, 1}, {1.0, 1.0});
  const ConicResult r = solve_conic(p, opt);
  CHECK(r.status == ConicStatus::Optimal);
  CHECK_NEAR(r.pobj, -1.5, 1e-7);
  CHECK_NEAR(r.dobj, -1.5, 1e-7);
  CHECK_NEAR(r.x[0], 0.0, 1e-6);
  CHECK_NEAR(r.x[1], 1.0, 1e-6);
  // With a cone: min 1/2 x^2 + t s.t. ||(x - 2)|| <= t (cone (t, x - 2)): x - 2 + ... on x < 2: 1/2 x^2 + 2 - x,
  // minimal at x = 1: 1.5. Variables (x, t); s = (t, x - 2) = h - G v with G = (0,-1), (-1,0), h = (0, -2).
  ConicProblem q = make({0, 1}, {}, {}, {{0, -1}, {-1, 0}}, {0, -2}, 0, {2});
  q.P = SparseMatrix::from_triplets(2, 2, {0}, {0}, {1.0});
  const ConicResult s = solve_conic(q, opt);
  CHECK(s.status == ConicStatus::Optimal);
  CHECK_NEAR(s.pobj, 1.5, 1e-7);
  CHECK_NEAR(s.x[0], 1.0, 1e-6);
  // Unbounded only along a direction P does not see: min 1/2 x^2 - y, y free -> dual infeasible.
  ConicProblem u = make({0, -1}, {}, {}, {{-1, 0}}, {10}, 1, {});
  u.P = SparseMatrix::from_triplets(2, 2, {0}, {0}, {1.0});
  CHECK(solve_conic(u, opt).status == ConicStatus::DualInfeasible);
}
