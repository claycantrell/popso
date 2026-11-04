#pragma once
// Local-MIP style local search (Lin, Zou, Cai: "New Local Search for Mixed Integer Programming",
// CP 2024 / AIJ 2025; the open-source Local-MIP 2.0 defaults). LP-free, works on rows normalised
// to a^T x <= b (equalities split in two):
//  - infeasible mode: "mixed tight moves" (one variable moves to the value that makes a violated
//    row tight / satisfied; integers rounded towards feasibility, continuous exact, clamped to the
//    bounds) on a BMS sample of violated rows, plus "breakthrough moves" (the same on the
//    objective row cost^T x <= best - delta once a solution is known); then tight moves on a
//    satisfied row (they use slack), a few "easy" moves (0 / bounds), and finally the best tight
//    move of one random violated row (after a weight update) whatever its score;
//  - scores are progress based (weight per row: +2w if it becomes satisfied, +w if a violated row
//    gets less violated, minus the same for the reverse) with an objective term once feasible, a
//    "breakthrough bonus" as tie-breaker, then the move's age; tabu on reversing a move;
//  - feasible mode: "lift moves": the objective variable that improves the objective most within
//    the range all rows allow (feasibility kept), maintained incrementally;
//  - monotone constraint weights (violated rows +1 at local minima; the objective row +1 when
//    feasible but not better than the best), restarts from the best point.
// Deterministic for a given seed and effort budget (effort = nonzeros touched, as in FJ).

#include <atomic>
#include <functional>
#include <vector>

#include "milp/search/mip_problem.hpp"

namespace milp {

struct LocalMipLimits {
  long long effort = 50000000;         // work units (nonzeros touched); deterministic budget
  // After this much effort without progress (a new best; before the first solution: fewer violated
  // rows than ever): stop, or, while no solution is known and stall_restart is set, restart from a
  // random point (the effort budget still applies).
  long long stall_effort = 100000000;
  bool stall_restart = false;
  double time_limit = kInf;            // seconds (safety cap)
  unsigned seed = 0;
  const std::atomic<bool>* stop = nullptr;  // external abort flag
  // Called every ~100k effort units with the effort so far (may block); false stops the search.
  std::function<bool(long long)> checkpoint;
  // Optional start point (e.g. an incumbent; clamped to the bounds, integers rounded). If it is
  // feasible it is the first solution (reported only if below `cutoff`).
  const std::vector<double>* start = nullptr;
  // Only points with cost^T x < cutoff (without offset) are reported; the search aims below it.
  double cutoff = kInf;
  long long restart_steps = 1000000;  // restart from the best point after this many steps without a new best
};

struct LocalMipStats {
  long long effort = 0, steps = 0, lift_moves = 0, restarts = 0;
  int solutions = 0;
  long long min_violated = -1;  // fewest violated (normalised) rows seen
  long long first_effort = -1;  // effort at the first reported solution
  std::vector<double> closest;  // without a solution: a point with the fewest violated rows
};

// Calls `on_solution(x, obj)` (obj = cost^T x, without offset) for each improving feasible point
// (max_violation <= 1e-6). Returns true if at least one was reported.
bool local_mip(const MipProblem& p, const LocalMipLimits& limits,
               const std::function<void(const std::vector<double>&, double)>& on_solution,
               LocalMipStats* stats = nullptr);

// Debug / research: MILP_LS_BENCH=<effort> runs FJ (restarted with new seeds until the budget is
// used) and local_mip on `p` with the same effort budget, single-threaded, and prints one line per
// method. Called on the presolved problem by the solver; returns false if the knob is not set.
bool local_search_bench(const MipProblem& p);

}  // namespace milp
