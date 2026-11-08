#pragma once
// Feasibility Jump (Luteberget & Sartor 2023): LP-free weighted violation
// local search. Each move "jumps" one variable to the value minimising the
// weighted violation of its rows; at local minima the weights of violated rows
// are bumped. After a feasible point is found an objective row
// (cost^T x <= best - delta) is added, so the search keeps improving it.

#include <atomic>
#include <functional>
#include <vector>

#include "milp/search/mip_problem.hpp"

namespace milp {

struct FjLimits {
  long long effort = 50000000;  // work units (nonzeros touched); deterministic budget
  long long stall_effort = 100000000;  // stop after this much effort without progress
  double time_limit = kInf;     // seconds (safety cap)
  unsigned seed = 0;
  const std::atomic<bool>* stop = nullptr;  // external abort flag
  // Called every ~1000 steps with the effort so far (may block); false stops FJ.
  std::function<bool(long long)> checkpoint;
  // > 0: after the first, the best point is reported at most every this many effort units and at the end (each
  // report is an incumbent check of the whole model: savsched1's 2335 improving points took 8.4 s of checks).
  long long report_effort = 0;
};

struct FjStats {
  long long effort = 0, moves = 0, restarts = 0;
  int solutions = 0;
  long long first_effort = -1;  // effort at the first solution
  long long min_violated = -1;  // fewest violated rows seen
  std::vector<double> closest;  // without a solution: a point with the fewest violated rows
};

// Calls `on_solution(x, obj)` (obj = cost^T x, without offset) for each
// improving feasible point found. Returns true if at least one was found.
bool feasibility_jump(const MipProblem& p, const FjLimits& limits,
                      const std::function<void(const std::vector<double>&, double)>& on_solution,
                      FjStats* stats = nullptr);

}  // namespace milp
