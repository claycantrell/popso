#pragma once
// Dual active-set method for the convex QPs of MIQP branch-and-bound nodes, and a safe
// Lagrangian lower bound for any primal / dual estimate.
//
//   min c^T x + 1/2 x^T Q x   s.t.  row_lb <= A x <= row_ub,  lb <= x <= ub
//
// The method is Goldfarb and Idnani's dual algorithm written on the KKT system of the working
// set: columns held at a bound are eliminated, and every iteration factors the dense reduced
// Hessian of the free columns (Cholesky) and the Schur complement of the active rows from
// scratch (no updates), so it pays off when few columns are free at the solution or the
// problem is small, which is the typical MIQP node. A working set that is optimal for the
// parent node is dual feasible for the child (whose bounds are tighter), so a child usually
// needs a few iterations. Q may be singular: the Hessian is regularized by a small proximal
// term eps_j (x_j - center_j)^2 / 2 around the start point, so the method solves a slightly
// perturbed QP; the bound reported is the safe bound of the true QP at the solution found
// (qp_safe_bound), so the perturbation (and any inaccuracy) can only weaken it.
//
// Safe bound: for any x and row multipliers y, convexity gives for every feasible point v
//   f(v) >= f(x) + g^T (v - x),  g = c + Q x,
// and g^T v = y^T A v + d^T v with d = g - A^T y is bounded below over the rows and the box
// (y_i > 0 uses row_lb_i, y_i < 0 row_ub_i; d_j > 0 uses lb_j, d_j < 0 ub_j). No optimality or
// feasibility of (x, y) is needed; reduced-cost fixing follows from the same terms.

#include <cstdint>
#include <memory>
#include <vector>

#include "milp/search/mip_problem.hpp"

namespace milp {

// Working set entries: column j at its lower bound 4j, at its upper bound 4j+1; row i at its
// lower side 4i+2, at its upper side 4i+3.
using WorkingSet = std::vector<int32_t>;

struct SafeBound {
  double bound = -kInf;
  std::vector<double> d;  // reduced costs c + Q x - A^T y (y cleaned: no multiplier on an infinite side)
};

// Q: lower triangle of the Hessian (with any convexity shift already added), qx = Q x (full
// symmetric product).
SafeBound qp_safe_bound(const MipProblem& p, const std::vector<double>& qx, const std::vector<double>& x,
                        const std::vector<double>& y, const std::vector<double>& lb, const std::vector<double>& ub);

struct ActiveSetResult {
  enum class Status { Optimal, Infeasible, Failed };
  Status status = Status::Failed;
  std::vector<double> x;  // Optimal: the point (bounded columns held at a bound exactly)
  std::vector<double> y;  // row multipliers, y_i >= 0 at an active lower side, <= 0 at an upper one
  SafeBound safe;         // safe bound of the true QP at (x, y)
  double objective = kInf;  // true objective of x (offset included, Q with shift)
  WorkingSet ws;
  uint64_t state = 0;  // id of the final factorization (kept for warm starts from ws), 0 if none
  int iterations = 0;
  long long flops = 0;
  const char* reason = "";  // why it failed
};

class ActiveSetQp {
 public:
  // diag_shift: added to Q's diagonal (may be empty). Builds a dense copy of Q (n x n).
  ActiveSetQp(const MipProblem& p, const std::vector<double>& diag_shift);
  ~ActiveSetQp();

  // start: working set to begin with (entries whose bound is infinite or that make the active
  // rows dependent are dropped, then constraints with a multiplier of the wrong sign); center:
  // the proximal center (size n; a parent's solution). flop_budget: the method gives up beyond.
  // start_state: the state id of the result whose ws is `start`; if that factorization is still
  // cached (the last few are), it is reused (the bound changes become borders) instead of
  // factoring the start working set again.
  ActiveSetResult solve(const std::vector<double>& lb, const std::vector<double>& ub, const WorkingSet& start,
                        const std::vector<double>& center, long long flop_budget, uint64_t start_state = 0);

  // Q x with Q including the shift (full symmetric product).
  void q_times(const std::vector<double>& x, std::vector<double>& out) const;
  // A working set guess from a primal point: bounds and row sides within tol (relative).
  WorkingSet guess(const std::vector<double>& x, const std::vector<double>& lb, const std::vector<double>& ub,
                   double tol = 1e-7) const;

 private:
  const MipProblem& p_;
  Index n_ = 0, m_ = 0;
  std::vector<double> Qd_;    // dense symmetric Q (+ shift), row-major n x n
  std::vector<double> eps_;   // proximal regularization per column
  SparseMatrix rows_;         // CSR copy of A
  std::vector<double> row_norm_;
  struct State;
  std::vector<std::unique_ptr<State>> cache_;  // most recent last
  uint64_t next_state_ = 1;
};

}  // namespace milp
