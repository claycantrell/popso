#pragma once
// Primal-dual interior-point method for second-order cone programs (v0.8):
//
//   min 1/2 x'Px + c'x  s.t.  A x = b,  G x + s = h,  s in K = R+^l x Q^{q_1} x ... x Q^{q_N}
//
// with P positive semidefinite (or empty), x free and Q^q = {(s0, s1) in R x R^{q-1} : s0 >= ||s1||}; the dual is
//
//   max -1/2 x'Px - b'y - h'z  s.t.  P x + A'y + G'z + c = 0,  z in K.
//
// The homogeneous self-dual embedding with Nesterov-Todd scaling and Mehrotra's predictor-corrector, as in ECOS
// (Domahidi, Chu, Boyd, "ECOS: An SOCP solver for embedded systems", ECC 2013) and CVXOPT's coneqp; with P the
// embedding's tau row carries x'Px / tau as in Clarabel (Goulart, Chen, "Clarabel: an interior-point solver for conic
// programs with quadratic objectives", 2024). Every iteration factors the regularized quasi-definite KKT matrix
//
//   [ P + delta I  A'       G'               ]
//   [ A         -delta I    0                ]
//   [ G          0         -W'W - delta I    ]
//
// (W'W diagonal on the orthant, one dense block per cone) by the sparse LDL' of lp/cholesky.hpp and solves with it
// three times (the tau direction, the affine and the combined direction), each with iterative refinement against
// the unregularized matrix. The embedding's tau and kappa give certificates: a primal infeasible problem ends with
// (y, z), h'z + b'y < 0 and A'y + G'z ~ 0; a dual infeasible (unbounded) one with (x, s), c'x < 0, A x ~ 0 and
// G x + s ~ 0. The data are equilibrated first (Ruiz: columns, rows of A, rows of G with one factor per cone).

#include <string>
#include <vector>

#include "milp/core/log.hpp"
#include "milp/core/sparse.hpp"
#include "milp/core/types.hpp"

namespace milp {

struct ConicProblem {
  Index n = 0;  // columns of x
  std::vector<double> c;
  SparseMatrix P;  // the objective's 1/2 x'Px: lower triangle with the diagonal, CSC (empty: none)
  SparseMatrix A;  // CSC, n columns, num_minor = rows (p)
  std::vector<double> b;
  SparseMatrix G;  // CSC, n columns, num_minor = rows (m = l + sum cones)
  std::vector<double> h;
  Index l = 0;               // the first l rows of G: the nonnegative orthant
  std::vector<Index> cones;  // then one second-order cone per entry, of that many rows (>= 1)
  double offset = 0.0;       // added to the objective values reported
};

enum class ConicStatus { Optimal, PrimalInfeasible, DualInfeasible, IterationLimit, TimeLimit, NumericalError };
const char* to_string(ConicStatus s);

struct ConicOptions {
  double feastol = 1e-8;  // primal and dual residuals relative to 1 + the data's norms
  double abstol = 1e-8;   // duality gap s'z, absolute
  double reltol = 1e-8;   // or relative to the objective
  int max_iterations = 100;
  double time_limit = kInf;  // seconds
  int threads = 1;
  int equilibration = 15;  // Ruiz passes (0: none)
  const Logger* log = nullptr;  // level 1: summary, 2: one line per iteration
};

struct ConicResult {
  ConicStatus status = ConicStatus::NumericalError;
  int iterations = 0;
  // Optimal (or at a limit): the point; infeasible: the certificate (y, z resp. x, s), normalized.
  std::vector<double> x, y, z, s;
  double pobj = 0.0, dobj = 0.0;  // 1/2 x'Px + c'x + offset, -1/2 x'Px - b'y - h'z + offset
  double pres = kInf, dres = kInf, gap = kInf;  // of the returned point, relative as in the stopping test
  // Optimal with residuals and gap only below 1e-6 (the best iterate once the method could not reach the options'
  // tolerances).
  bool reduced_accuracy = false;
  double time = 0.0;
};

ConicResult solve_conic(const ConicProblem& prob, const ConicOptions& opt);

}  // namespace milp
