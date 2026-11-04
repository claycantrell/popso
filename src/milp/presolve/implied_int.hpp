#pragma once
// Implied-integer detection (continuous columns whose integrality comes for free).
//
// Rules, iterated to a fixed point (columns are added one at a time, each justified by integer
// columns, integrally fixed columns and columns found before it):
//   - primal: an equation sum a_k x_k = b in which x_j is the only column that is neither integer nor
//     implied integer, with every a_k / a_j and b / a_j integral: x_j is integral whenever the others are.
//   - dual (optional): x_j has integral (or infinite) bounds and, in every row of x_j, the other columns
//     are (implied) integer with integral a_k / a_j and the finite sides integral multiples of a_j. For
//     fixed values of the others, x_j's feasible set is then an interval with integral end points, so any
//     feasible point can move x_j to an integral end point without worsening the (linear) objective.
// Because a dual column's rows contain only columns found before it, the rules compose: some optimal
// solution has every implied integer integral (SCIP's / HiGHS's implied integer detection).

#include <vector>

#include "milp/search/mip_problem.hpp"

namespace milp {

struct ImpliedIntStats {
  Index primal = 0, dual = 0;  // newly detected columns by rule
};

// Flags per column (1 = implied integer; integer columns are 0). Columns already flagged in
// p.implied_int count as implied integers. bound_tol: a finite bound within this distance of an
// integer counts as integral for the dual rule (the caller snaps it). The dual rule is skipped
// when `dual` is false, for quadratic objectives and with logical implications.
std::vector<char> detect_implied_integers(const MipProblem& p, bool dual = true, double bound_tol = 1e-9,
                                          ImpliedIntStats* stats = nullptr);

// Detects implied integers of p (keeping its existing flags), rounds their finite bounds to integers
// (snapping bounds within 1e-6 of an integer; a flag whose rounded bounds would cross is dropped)
// and stores the flags in p.implied_int. Returns the number of flagged columns.
Index mark_implied_integers(MipProblem& p, bool dual = true, ImpliedIntStats* stats = nullptr);

}  // namespace milp
