#pragma once
// Branch-and-bound for mixed-integer convex quadratic programs (convex objective
// c^T x + 1/2 x^T Q x, linear constraints, integer columns).
//
// The root relaxation is solved by the barrier (qp.hpp); the nodes by the dual active-set method
// of qp_active.hpp, warm-started from the parent's working set (and its cached factorization),
// with the barrier as the fallback. Node bounds are safe Lagrangian bounds of the relaxation
// (valid for any primal point and multipliers). Single-threaded and deterministic:
//  - node selection: best bound with plunging (the child with the smaller pseudocost estimate,
//    or the rounding direction, is processed next, its sibling is queued; a plunge ends at a
//    pruned, infeasible or integral node);
//  - branching: reliability branching on the QP objective (product score; candidates with
//    fewer than 4 pseudocost observations are strong-branched by the active-set method, a
//    child that is infeasible or cut off fixes the column to the other side);
//  - domain propagation of the linear rows (search/propagation.hpp) at every node and
//    reduced-cost tightening of integer bounds from the safe bound (at the nodes, and globally
//    from the root whenever the incumbent improves); no cutting planes;
//  - primal heuristics: simple rounding at every node of a pure-integer problem, rounding with
//    the continuous columns re-optimized (integers fixed, QP) at the root, a QP dive (fix the
//    least fractional variable and the integral ones, re-solve, one backtrack per level) at
//    the root and periodically, and Feasibility Jump on the MILP with the objective
//    linearized at the root relaxation's solution (its integer assignments are completed by a
//    QP over the continuous columns). Every incumbent is evaluated with the true objective.
//  - a node whose relaxation looks integral but whose rounded point cannot be completed to a
//    feasible one is split three ways on an integer column (x <= v-1, x = v, x >= v+1).

#include <functional>
#include <vector>

#include "milp/core/log.hpp"
#include "milp/core/timer.hpp"
#include "milp/search/mip_problem.hpp"
#include "milp/solver/solver.hpp"

namespace milp {

struct MiqpOptions {
  int threads = 1;                  // for the factorizations of large relaxations
  std::vector<double> diag_shift;   // added to Q's diagonal in the relaxations (convexity repair)
  std::vector<double> start;        // starting point (one value per column, NaN = not given) or empty
  // Each new incumbent: objective in the user's sense, what found it, time (SolverParams::trace_incumbents).
  std::function<void(double objective, const char* source, double time)> on_trace;
};

struct MiqpResult {
  Status status = Status::Unknown;
  bool has_solution = false;
  std::vector<double> x;
  double objective = kInf;  // of p (minimization, offset included)
  double bound = -kInf;
  double root_bound = -kInf;
  long long nodes = 0;
  long long barrier_iterations = 0;
  long long lp_iterations = 0;
};

// p.Q must be convex (see check_convexity). params: time and node limits, gaps, verbosity
// (through log), memory limit.
MiqpResult miqp_solve(const MipProblem& p, const SolverParams& params, const Logger& log, const Timer& timer,
                      const MiqpOptions& opt);

}  // namespace milp
