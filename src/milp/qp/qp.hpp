#pragma once
// Convex quadratic programming: convexity test of the Hessian and the QP solve
//
//   min c^T x + 1/2 x^T Q x + offset   s.t.  row_lb <= A x <= row_ub,  lb <= x <= ub
//
// by the barrier (interior-point) method, lp/barrier.hpp. There is no basis: the result is the
// primal point (and the objective); infeasibility and unboundedness reported by the barrier
// are confirmed by the dual simplex on auxiliary LPs (feasibility of the rows and bounds, and
// a recession direction d with Q d = 0, c^T d < 0), so they are never reported on a hunch.

#include <atomic>
#include <string>
#include <vector>

#include "milp/core/log.hpp"
#include "milp/core/sparse.hpp"
#include "milp/core/types.hpp"
#include "milp/search/mip_problem.hpp"

namespace milp {

struct ConvexityResult {
  bool convex = true;
  // Q was positive semidefinite only up to a small tolerance: diag_shift (per column, may be
  // empty) is added to Q's diagonal for the solve to make the relaxation convex. The objective
  // reported is always computed with the original Q.
  bool perturbed = false;
  std::vector<double> diag_shift;
  std::string message;  // why the objective is not convex (or what was perturbed)
};

// Q: lower triangle (row >= column) of the Hessian of a minimization, CSC n x n. The test
// scales Q to a unit diagonal (D^-1/2 Q D^-1/2; a zero diagonal entry with a nonzero
// off-diagonal one in its column, or a negative diagonal entry, is nonconvex at once) and
// factors it plus 1e-10 I by the sparse Cholesky: a negative pivot means a negative eigenvalue.
// If only the shift 1e-6 I succeeds, the objective is accepted as convex with that relative
// perturbation of the diagonal (perturbed = true).
ConvexityResult check_convexity(const SparseMatrix& Q);

enum class QpStatus { Optimal, Infeasible, Unbounded, InfeasibleOrUnbounded, TimeLimit, NumericalError };
const char* to_string(QpStatus s);

struct QpOptions {
  double tol = 1e-8;  // barrier tolerance (relative primal/dual infeasibility and gap)
  double time_limit = kInf;
  int threads = 1;
  const Logger* log = nullptr;  // barrier log (level 1 summary, 2 iterations)
  const std::atomic<bool>* stop = nullptr;
  bool confirm = true;  // confirm infeasibility / unboundedness by simplex LPs
  int max_iterations = 300;
};

struct QpSolution {
  QpStatus status = QpStatus::NumericalError;
  std::vector<double> x;    // primal point (Optimal), clipped to the bounds
  // Row duals of the barrier's final iterate (Optimal): the reduced costs are
  // d = c + Q x - A^T y (qp_reduced_costs), y_i >= 0 at an active lower row bound, <= 0 at an upper one.
  std::vector<double> y;
  double objective = kInf;  // offset + c^T x + 1/2 x^T Q x of x
  double dual_bound = -kInf;  // barrier dual objective estimate (offset included)
  double primal_infeas = kInf;  // max row / bound violation of x
  int iterations = 0;          // barrier iterations (all attempts)
  bool approximate = false;    // optimal only to about 1e-6
};

// d = c + Q x - A^T y of p (Q from its lower triangle).
std::vector<double> qp_reduced_costs(const MipProblem& p, const std::vector<double>& x, const std::vector<double>& y);

// Solves the continuous relaxation of p (integrality ignored) with column bounds lb, ub and
// the Hessian p.Q plus the optional diagonal shift (size p.num_cols() or empty).
QpSolution solve_qp(const MipProblem& p, const std::vector<double>& lb, const std::vector<double>& ub,
                    const QpOptions& opt, const std::vector<double>& diag_shift = {});

}  // namespace milp
