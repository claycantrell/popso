#pragma once
// Clique merging: replaces set-packing rows (sum of binaries <= 1) by fewer,
// larger clique rows covering them, e.g. the pairwise conflicts
// x_a + x_b <= 1 of a pigeonhole model become one row per hole. Only rows
// change, so no postsolve is needed.

#include "milp/search/mip_problem.hpp"

namespace milp {

// Returns the number of set-packing rows removed (dominated by new cliques).
Index merge_cliques(MipProblem& p, Index max_new_rows = 10000);

}  // namespace milp
