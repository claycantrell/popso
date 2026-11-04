#pragma once
// Dual values of LPs and fixed-MIP duals, recovered in the original (unpresolved) space.
//
// Presolve only maps primal solutions back (presolve/presolve.hpp); several of its reductions
// (dual fixing, dominated columns, parallel columns, substitutions) have no simple dual postsolve.
// So the duals are not carried through presolve. Instead, the optimal primal point x (postsolved
// to the original columns) is turned into an optimal basis of the original LP: a crossover from x
// (lp/crossover.hpp, primal push only) gives a primal feasible basis near x, and the simplex
// finishes it (primal simplex first, then the dual simplex). Duals of an optimal basis are optimal,
// and any optimal primal-dual pair satisfies complementary slackness, so the duals are
// complementary to x itself (up to the tolerances), not only to the basic solution.
//
// For a MILP the "fixed-MIP" duals are those of the LP with every integer column fixed at its
// (rounded) incumbent value.

#include <cstdint>
#include <vector>

#include "milp/lp/simplex.hpp"
#include "milp/search/mip_problem.hpp"

namespace milp {

struct LpDualResult {
  bool ok = false;
  std::vector<double> row_dual;      // y (minimization of p: y_i >= 0 at an active row lower bound)
  std::vector<double> reduced_cost;  // d = c - A^T y
  // Basis status of the optimal basis: columns and rows (a row is "at lower" when its activity is
  // at row_lb; VarStatus::AtZero is a free nonbasic column or row).
  std::vector<VarStatus> col_status, row_status;
  long long iterations = 0;  // simplex iterations of the cleanup
  double objective = 0.0;    // c^T x of the basic solution, plus p.offset
};

// p is solved as an LP (integrality ignored) with column bounds col_lb / col_ub (p's own if
// empty); x: an optimal (or, for fixed-MIP duals, feasible) point of that LP.
LpDualResult recover_lp_duals(const MipProblem& p, const std::vector<double>& col_lb, const std::vector<double>& col_ub,
                              const std::vector<double>& x, double time_limit);

// Consistency of recovered duals: the objective of the basis they come from (LpDualResult::
// objective) against the solution's, relative tolerance 1e-6. For an LP both are the optimum
// (Consistent); for fixed-MIP duals the incumbent is feasible for the fixed LP, so the basis may
// also be better (BetterThanSolution: the incumbent's continuous part is not optimal for the fixed
// LP, the duals belong to the fixed LP's optimum). Anything else is a failed recovery (Mismatch:
// the duals are dropped).
enum class DualObjectiveCheck { Consistent, BetterThanSolution, Mismatch };
DualObjectiveCheck check_dual_objective(double basis_objective, double solution_objective, bool fixed_mip);

}  // namespace milp
