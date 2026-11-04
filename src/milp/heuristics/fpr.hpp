#pragma once
// Fix-propagate-repair (Salvagnin, Roberti, Fischetti, "A fix-propagate-repair heuristic for mixed
// integer programming", Math. Prog. Comp. 17 (2025) 111-139).
//
// A limited depth-first search over the integer variables: they are fixed one at a time in a static
// order to a value chosen by a rule, each fixing is propagated (DomainPropagator), and when a row
// becomes violated a WalkSAT-like repair walk shifts already fixed variables (and the domains of
// continuous ones) to restore feasibility of the partial assignment before the search backtracks.
// At the bottom of a dive the continuous variables get their values from an LP with the integers
// fixed. Single strategies are not robust, so a portfolio of variable-order / value-rule /
// search-mode combinations is run until one finds a feasible point. LP-free unless an LP point is
// given (it enables the LP-guided value rule). Deterministic for a given seed and work budget.

#include <atomic>
#include <functional>
#include <vector>

#include "milp/search/mip_problem.hpp"

namespace milp {

// Search parametrizations of the paper (Sect. 3): which of propagation, repair and backtracking
// on infeasible nodes are on.
enum class FprMode : uint8_t {
  Dfs,       // propagate, no repair, backtrack (plain fix-and-propagate)
  DfsRep,    // propagate, repair, backtrack
  Dive,      // no propagation, repair, no backtrack
  DiveProp,  // propagate, repair, no backtrack
};

// Static variable orders (Table 1); continuous variables are never branched on.
enum class FprOrder : uint8_t {
  LeftRight,    // formulation order
  Type,         // binaries, then general integers (formulation order within each)
  Random,       // as Type, shuffled within each type
  Locks,        // as Type, by non-increasing max(up locks, down locks) within each type
  TypeCliques,  // as Type, binaries grouped by a greedy clique cover
};

// Value rules (Table 2); the value is always clipped to the current domain.
enum class FprValue : uint8_t {
  Up,        // upper bound
  Random,    // lower or upper bound at random
  GoodObj,   // the bound the objective prefers (lower bound if the cost is 0)
  BadObj,    // the other bound (lower bound if the cost is 0)
  LooseDyn,  // the direction with fewer locks among the rows not yet redundant
  Lp,        // round the reference LP value up with probability equal to its fractional part
};

struct FprStrategy {
  FprMode mode;
  FprOrder order;
  FprValue value;
  const char* name;
};

// The LP-free portfolio of the paper (Sect. 6.3) and the LP-guided one (needs FprLimits::lp_x).
const std::vector<FprStrategy>& fpr_lp_free_portfolio();
const std::vector<FprStrategy>& fpr_lp_portfolio();

struct FprLimits {
  long long effort = 1000000000;  // work units (nonzeros touched) for the whole call; deterministic budget
  long long run_effort = 0;       // per strategy run (0: 100 * nnz + 100000, as in the paper)
  double time_limit = kInf;       // seconds (safety cap; makes the result timing dependent if it binds)
  unsigned seed = 0;
  int rounds = 0;                 // passes over the portfolio (randomized strategies get fresh seeds); 0: until the budget is spent
  const std::atomic<bool>* stop = nullptr;
  // Reference LP point (columns of p): runs the LP-guided strategies after the LP-free ones.
  const std::vector<double>* lp_x = nullptr;
  std::vector<FprStrategy> strategies;  // empty: the portfolio(s) above
  bool one_opt = true;                  // improve the point found with 1-opt shifts
};

struct FprStats {
  long long effort = 0;
  int runs = 0;  // strategy runs started
  long long nodes = 0, repairs = 0, repairs_ok = 0, leaf_lps = 0, leaf_lps_ok = 0;
  const char* strategy = nullptr;  // the strategy that found the point
  bool infeasible = false;         // root propagation proved the problem infeasible
};

// Runs the portfolio until a strategy finds a feasible point (max violation <= 1e-6), written to
// `sol`. Returns false if none was found within the limits.
bool fix_propagate_repair(const MipProblem& p, const FprLimits& limits, std::vector<double>& sol,
                          FprStats* stats = nullptr);

class Logger;
// Standalone comparison of FJ, local MIP and FPR on the same work budget (MILP_FPR_BENCH):
// prints one "HEURBENCH" line per heuristic with the time and the effort to the first feasible
// point, the first and the best objective (user's sense, incl. the offset), and "fpr-lmip" (local
// MIP started from FPR's point with the rest of the budget). Knobs: MILP_FPR_BENCH_NOFJ /
// _NOLMIP skip a heuristic, _LP adds the LP-guided portfolio after the root LP, _FJ FPR guided by
// FJ's closest point, _EACH every LP-free strategy alone.
void heuristic_bench(const MipProblem& p, long long effort, unsigned seed, const Logger& log);

}  // namespace milp
