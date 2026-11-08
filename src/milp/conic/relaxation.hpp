#pragma once
// The continuous relaxation of a MipProblem with convex quadratic rows (MipProblem::quad_rows) as a second-order
// cone program (conic/socp.hpp), over the same columns (integrality dropped):
// - finite bounds and one-sided or ranged linear rows on the nonnegative orthant, equalities (and fixed columns)
//   in A;
// - a convex quadratic row a'x + x'Sx <= hi (S PSD; a finite lower side is negated) as the rotated cone
//   ||F x||^2 <= hi - a'x, written as the cone ((1 + u) / 2, (u - 1) / 2, F x) with u = hi - a'x, where
//   S = F'F by a pivoted Cholesky factorization of S on the row's columns (a diagonal S directly);
// - a cone row (QuadRow::cone) sqrt(gamma + sum (a_m x + b_m)^2) <= h x + h0 as the cone
//   (h x + h0, sqrt(gamma), a_1 x + b_1, ...);
// - the objective's quadratic part (MipProblem::Q, or the epigraph row of an outer approximation) as P; the epigraph
//   column is then fixed at 0 (the caller sets it to 1/2 x'Qx at the solution).

#include <string>

#include "milp/conic/socp.hpp"
#include "milp/search/mip_problem.hpp"

namespace milp {

// Where each row of the conic problem came from, to map its duals back: A rows and orthant G rows are sign times a
// MipProblem row (kind 0) or a column bound (kind 1), or neither (kind 2: the epigraph column fixed at 0); each cone
// is a quadratic row (its index in quad_rows, its direction: +1 for an upper side) of the rotated form (norm = false:
// the multiplier is (z0 + z1) / 2) or the norm form of a cone row (norm = true: z0 / (2 (h x + h0))).
struct ConicMap {
  std::vector<int8_t> a_kind, g_kind;
  std::vector<Index> a_index, g_index;
  std::vector<double> a_sign, g_sign;
  std::vector<Index> cone_quad;
  std::vector<double> cone_dir;
  std::vector<char> cone_norm;
};

// False, with the reason in *why, if a row's factorization finds it not positive semidefinite.
bool conic_relaxation(const MipProblem& p, ConicProblem& out, std::string* why = nullptr, ConicMap* map = nullptr);

}  // namespace milp
