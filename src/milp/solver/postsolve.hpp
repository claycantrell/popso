#pragma once
// The search's answer on the presolved problem, mapped back to the original columns.

#include <vector>

#include "milp/core/log.hpp"
#include "milp/core/timer.hpp"
#include "milp/presolve/presolve.hpp"
#include "milp/search/branch_and_bound.hpp"
#include "milp/search/mip_problem.hpp"
#include "milp/solver/solver.hpp"

namespace milp {

// Re-optimizes the continuous columns of x by an LP with the integer columns fixed at their rounded
// values (unperturbed, primal feasibility tolerance primal_tol). Returns false, leaving x as it is, if
// that LP fails or its point violates p by more than 1e-6.
bool polish_point(const MipProblem& p, const Timer& timer, double time_limit, std::vector<double>& x,
                  double primal_tol = 1e-6);

// Postsolves mr (the search's answer on pre.reduced(): x in its columns, objective and bound with its
// offset) into p's columns and sets mr->objective to the objective of the point in p.
//
// A point that violates p by more than 1e-6 is repaired by polish_point. A point whose objective came
// out worse than the search's gets its continuous columns re-optimized at a tight tolerance: postsolve
// keeps the objective up to round-off except where it clamps a substituted column into its bounds
// (see the implementation).
//
// Returns false (and clears the solution) if the point cannot be repaired: p should then be solved
// without presolve.
bool postsolve_result(const MipProblem& p, const Presolver& pre, const SolverParams& params, const Timer& timer,
                      const Logger& log, MipResult* mr);

}  // namespace milp
