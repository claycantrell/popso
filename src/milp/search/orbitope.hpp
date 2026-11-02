#pragma once
// Packing/partitioning orbitopes (Kaibel & Pfetsch) and orbitopal fixing.
//
// A full orbitope is an r x p matrix of binary columns on which the symmetry
// group acts as the symmetric group S_p permuting its columns (e.g. identical
// trucks, machines or bins: x[task][truck]). Detection follows HiGHS: a
// component of the generators found by detect_symmetry whose generators are
// all involutions with r 2-cycles each, and each either swaps two existing
// matrix columns or maps one existing column onto a new one.
//
// If every matrix row is a set packing row (at most one 1, e.g. a set
// packing/partitioning constraint of the model), the representatives with
// lexicographically decreasing columns are exactly the matrices whose
// columns' leading ones alpha_0 < alpha_1 < ... (rows of their first 1) are
// strictly increasing, all-zero columns last. Restricting the search to these
// representatives keeps an optimal solution. Orbitopal fixing derives from the
// local 0/1 fixings:
//   lower bounds on alpha_c (first row not fixed to 0, after alpha_{c-1}):
//     x_ic = 0 above it, the whole column if there is none;
//   upper bounds on alpha_c (first row fixed to 1, or alpha_{c+1} - 1):
//     infeasible if below the lower bound, x_ic = 1 if both coincide.
// The search must then not use other symmetry reductions on these columns
// (orbital branching uses the remaining generators only).

#include <vector>

#include "milp/search/mip_problem.hpp"
#include "milp/search/symmetry.hpp"

namespace milp {

struct Orbitope {
  Index rows = 0, cols = 0;
  std::vector<Index> var;  // var[i * cols + c]: column of the model at matrix entry (i, c)
  Index at(Index i, Index c) const { return var[static_cast<size_t>(i) * static_cast<size_t>(cols) + c]; }
  // swaps[c]: the verified symmetry exchanging matrix columns c and c + 1 (moved columns and images).
  std::vector<std::vector<std::pair<Index, Index>>> swaps;
};

// Maps a point onto its representative (lexicographically decreasing matrix columns) by
// exchanging adjacent matrix columns with the verified swaps (debugging: MILP_DEBUG_SOL).
void orbitope_canonicalize(const Orbitope& o, std::vector<double>& x);

// Packing orbitopes among the generators of sym. used_gen[g] is set for every
// generator of a component that became an orbitope.
std::vector<Orbitope> detect_packing_orbitopes(const MipProblem& p, const Symmetry& sym,
                                               std::vector<char>* used_gen);

// Orbitopal fixing on local bounds [lb, ub]: appends each tightened column to
// `tightened` (and its new bound to `fixed_one` as 1 if fixed to one, 0 if to
// zero). Returns false if no representative lies within the bounds.
bool orbitope_propagate(const Orbitope& o, std::vector<double>& lb, std::vector<double>& ub,
                        std::vector<Index>& tightened, std::vector<char>& fixed_one);

}  // namespace milp
