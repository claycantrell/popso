#pragma once
// Large-neighbourhood heuristics that solve a restricted copy of the problem
// with a small branch-and-bound (no recursive sub-MIPs).

#include <atomic>
#include <vector>

#include "milp/search/mip_problem.hpp"

namespace milp {

// Outcome of the sub-MIPs solved by one heuristic call (accumulated over retries).
struct SubMipStats {
  int solves = 0;              // sub-MIPs actually solved (0: neighbourhood rejected)
  bool exhausted = false;      // last sub-MIP finished (optimal or infeasible): neighbourhood too small
  long long lp_iterations = 0, nodes = 0;
  Index fixed = 0, ints = 0;  // integer variables fixed / total in the last sub-MIP
};

struct SubMipLimits {
  long long node_limit = 500;
  double time_limit = 5.0;  // seconds
  double min_fixed_fraction = 0.3;  // skip if fewer integer variables are fixed
  SubMipStats* stats = nullptr;     // optional output
  std::atomic<bool>* stop = nullptr;           // optional: abort when set (result: whatever was found)
  std::atomic<long long>* progress = nullptr;  // optional: LP iterations so far are added here
  long long iteration_limit = -1;  // LP iterations per sub-MIP solve (-1: none); deterministic, unlike time
  int mem_limit_mb = 256;          // search memory cap of the sub-MIP
};

// Solves p restricted to column bounds [lb, ub]. Returns true if a solution
// with internal objective (incl. offset) below `cutoff` was found.
bool solve_sub_mip(const MipProblem& p, const std::vector<double>& lb, const std::vector<double>& ub, double cutoff,
                   const SubMipLimits& limits, std::vector<double>& sol);

// RENS: fix integers that are integral in the LP solution, restrict the rest to floor/ceil.
// With `redcost` given, an infeasible sub-MIP is retried fixing fewer variables
// (keeping those with the largest |reduced cost| fixed).
bool rens(const MipProblem& p, const double* lp_x, double cutoff, const SubMipLimits& limits,
          std::vector<double>& sol, const std::vector<double>* redcost = nullptr);

// RINS: fix integers on which the LP solution and the incumbent agree. If that fixes less
// than `target_fixed` of the integers, more are fixed at their incumbent values, those
// whose LP value is closest to the incumbent first.
bool rins(const MipProblem& p, const double* lp_x, const std::vector<double>& incumbent, double cutoff,
          const SubMipLimits& limits, std::vector<double>& sol, double target_fixed = 0.0);

// Repair: fix the integers of an infeasible point x except those in its violated
// rows and in rows sharing variables with them (at most about half of the integers stay free).
bool repair(const MipProblem& p, const std::vector<double>& x, double cutoff, const SubMipLimits& limits,
            std::vector<double>& sol);

// Mutation: fix a random `fix_fraction` of the integer variables (that are not
// fixed already) to their incumbent values.
bool mutation(const MipProblem& p, const std::vector<double>& incumbent, double fix_fraction, unsigned seed,
              double cutoff, const SubMipLimits& limits, std::vector<double>& sol);

// Crossover: fix the integer variables on which all the given feasible points agree (their
// common value); the others keep their bounds. Needs at least two points.
bool crossover(const MipProblem& p, const std::vector<const std::vector<double>*>& points, double cutoff,
               const SubMipLimits& limits, std::vector<double>& sol);

}  // namespace milp
