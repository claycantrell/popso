#pragma once
// Public entry point of the solver.

#include <atomic>
#include <cstdint>
#include <functional>
#include <string>
#include <vector>

#include "milp/core/model.hpp"
#include "milp/core/types.hpp"

namespace milp {

enum class Status {
  Unknown,
  Optimal,
  Infeasible,
  Unbounded,
  InfeasibleOrUnbounded,
  TimeLimit,
  NodeLimit,
  MemoryLimit,
  IterationLimit,
  Interrupted,
  NumericalError,
  NonConvex,  // the quadratic objective is not convex (minimization) / concave (maximization)
};

const char* to_string(Status s);

// Method for the LP of a pure LP model and for the root LP of a MIP. Auto picks the barrier for
// large LPs where it was measured to win (see branch_and_bound.cpp, use_barrier_root).
enum class LpMethod { Auto = 0, Simplex = 1, Barrier = 2 };
const char* to_string(LpMethod m);

struct SolverParams {
  double time_limit = kInf;   // seconds
  int threads = 0;            // 0 = hardware concurrency
  int verbosity = 1;          // 0 silent, 1 summary, 2 progress, 3 debug
  double mip_rel_gap = 1e-4;  // relative gap termination
  double mip_abs_gap = 1e-10;
  long long node_limit = -1;  // -1 = unlimited
  long long mem_limit_mb = 4096;  // approximate cap on search-tree memory
  unsigned seed = 0;
  int racing = 1;             // concurrent differently-seeded searches sharing incumbents (opportunistic mode)
  bool deterministic = false;  // reproducible parallel search (same seed + threads -> same result); slower
  bool presolve = true;
  // Symmetry handling: detection, symmetry-breaking rows, orbitopal fixing, orbital branching.
  bool symmetry = true;
  LpMethod lp_method = LpMethod::Auto;
  // Logical constraints (semi-continuous, SOS, indicators): largest big-M coefficient of their
  // reformulation; beyond it (or without finite bounds) they are enforced by branching. 0: always branching.
  double logic_big_m = 1e6;
  // Dual values (SolveResult::row_dual etc.): 0 none, 1 for LPs, 2 also the fixed-MIP duals of a
  // MILP (duals of an optimal basis of the LP with the integer columns fixed at the incumbent; after
  // an early stop they need not be complementary to the incumbent's continuous part, which is left
  // as it is). Duals failing the objective consistency check are dropped. Continuous QPs get
  // their duals in any case (they come with the barrier's solution). See solver/duals.hpp.
  int duals = 0;
  // Interrupt flag (core/interrupt.hpp): once set (from any thread or a signal handler), the solve
  // stops soon with Status::Interrupted and the best solution found so far. Neither reset nor set by
  // solve() (a stop requested by on_incumbent uses an internal flag), so params can be reused.
  std::atomic<bool>* interrupt = nullptr;
  // Called for each new incumbent, in the model's columns: its objective and the current dual
  // bound (both in the model's sense). Calls are serialized and the objectives strictly improve;
  // in deterministic mode they come at deterministic points of the search. Returning true asks the
  // solve to stop (as an interrupt: Status::Interrupted). Not for quadratic objectives.
  std::function<bool(const std::vector<double>& x, double objective, double bound)> on_incumbent;
  Tolerances tol;
};

struct SolveStats {
  double time = 0.0;
  long long nodes = 0;
  long long lp_iterations = 0;
  double root_bound = -kInf;
  long long barrier_iterations = 0;  // interior-point iterations (QP and MIQP relaxations)
};

struct SolveResult {
  Status status = Status::Unknown;
  bool has_solution = false;
  double objective = kInf;   // in the model's original sense
  double bound = -kInf;      // best dual bound, original sense
  std::vector<double> x;     // primal solution (if has_solution)
  // Duals (continuous QPs and, with SolverParams::duals, LPs solved to optimality and fixed-MIP
  // duals of MILPs with a solution; empty otherwise): row duals y and reduced costs
  // d = grad f(x) - A^T y with f the model's objective (c + Q x). Minimization: y_i >= 0 at an
  // active row lower bound, <= 0 at an upper one, d_j >= 0 at a column lower bound, <= 0 at an
  // upper one, 0 in between; maximization: the opposite signs.
  std::vector<double> row_dual;
  std::vector<double> reduced_cost;
  bool has_duals = false;  // row_dual and reduced_cost are set (the model may have no rows)
  // LPs and fixed-MIP duals (SolverParams::duals) also give the optimal basis they come from:
  // 0 basic, 1 at lower bound, 2 at upper bound, 3 nonbasic free (zero); a row is at its lower
  // bound when its activity is at row_lb. Empty otherwise (continuous QPs have no basis).
  std::vector<int8_t> col_basis, row_basis;
  bool fixed_mip_duals = false;  // the duals are those of the MILP with its integers fixed
  SolveStats stats;
};

// Models with a quadratic objective (Model::q_*) are presolved by qp/qp_presolve.hpp and solved as
// convex QPs by the barrier method (no basis: x, the objective and the duals; see qp/qp.hpp) or,
// with integer columns, by the MIQP branch-and-bound of qp/miqp.hpp. A nonconvex objective gives
// Status::NonConvex.
SolveResult solve(const Model& model, const SolverParams& params);

// Solves with a starting solution: one value per model column, NaN = not given
// (see start.hpp). An empty vector means no start.
SolveResult solve(const Model& model, const SolverParams& params, const std::vector<double>& start);

}  // namespace milp
